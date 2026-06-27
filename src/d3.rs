//! D3Q19 lattice-Boltzmann — the 3D extension that yields a TRUE vehicle drag
//! coefficient (force normalised by frontal area), not just the 2D sectional
//! ranking. Same BGK scheme as the 2D solver, one dimension up. A car-ish body
//! is voxelised; drag is measured by momentum exchange and divided by the
//! projected frontal cell count.

/// D3Q19 velocities.
const E: [(i32, i32, i32); 19] = [
    (0, 0, 0),
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
    (1, 1, 0), (-1, -1, 0), (1, -1, 0), (-1, 1, 0),
    (1, 0, 1), (-1, 0, -1), (1, 0, -1), (-1, 0, 1),
    (0, 1, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1),
];
const W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
];
const OPP: [usize; 19] =
    [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17];

/// A voxelised 3D body in the wind-tunnel. z is height (ground at z=0), flow +x.
pub struct Body3 {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub solid: Vec<bool>,
    /// Projected frontal area (count of y-z columns containing any solid cell).
    pub frontal_cells: usize,
}

impl Body3 {
    fn frontal(nx: usize, ny: usize, nz: usize, solid: &[bool]) -> usize {
        let mut c = 0;
        for z in 0..nz {
            for y in 0..ny {
                let mut hit = false;
                for x in 0..nx {
                    if solid[(z * ny + y) * nx + x] { hit = true; break; }
                }
                if hit { c += 1; }
            }
        }
        c
    }

    /// A squareback box car: full-height block, leading edge at `x0`.
    pub fn box_car(nx: usize, ny: usize, nz: usize, x0: usize, len: usize, w: usize, h: usize) -> Body3 {
        let mut solid = vec![false; nx * ny * nz];
        let y0 = (ny - w) / 2;
        for x in x0..(x0 + len).min(nx) {
            for y in y0..(y0 + w).min(ny) {
                for z in 1..(1 + h).min(nz) {
                    solid[(z * ny + y) * nx + x] = true;
                }
            }
        }
        let frontal_cells = Self::frontal(nx, ny, nz, &solid);
        Body3 { nx, ny, nz, solid, frontal_cells }
    }

    /// Same frontal face, but the roof tapers down over the rear half
    /// (fastback) — lower drag at equal frontal area.
    pub fn fastback_car(nx: usize, ny: usize, nz: usize, x0: usize, len: usize, w: usize, h: usize) -> Body3 {
        let mut solid = vec![false; nx * ny * nz];
        let y0 = (ny - w) / 2;
        let half = len / 2;
        for k in 0..len {
            let x = x0 + k;
            if x >= nx { break; }
            // full height for the front half, tapering to ~40% at the tail
            let hz = if k < half {
                h
            } else {
                let frac = 1.0 - 0.6 * ((k - half) as f64) / (half.max(1) as f64);
                ((h as f64) * frac).round() as usize
            };
            for y in y0..(y0 + w).min(ny) {
                for z in 1..(1 + hz).min(nz) {
                    solid[(z * ny + y) * nx + x] = true;
                }
            }
        }
        let frontal_cells = Self::frontal(nx, ny, nz, &solid);
        Body3 { nx, ny, nz, solid, frontal_cells }
    }

    #[inline]
    fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        self.solid[(z * self.ny + y) * self.nx + x]
    }
}

pub struct Lbm3 {
    nx: usize,
    ny: usize,
    nz: usize,
    tau: f64,
    u0: f64,
    f: Vec<f64>,
    ftmp: Vec<f64>,
    body: Body3,
}

#[inline]
fn feq(i: usize, rho: f64, u: (f64, f64, f64)) -> f64 {
    let (ex, ey, ez) = E[i];
    let eu = ex as f64 * u.0 + ey as f64 * u.1 + ez as f64 * u.2;
    let usq = u.0 * u.0 + u.1 * u.1 + u.2 * u.2;
    W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq)
}

impl Lbm3 {
    pub fn new(body: Body3, re: f64) -> Lbm3 {
        let (nx, ny, nz) = (body.nx, body.ny, body.nz);
        let u0 = 0.05;
        // reference length = sqrt(frontal area) for the Reynolds/viscosity set
        let d = (body.frontal_cells as f64).sqrt().max(1.0);
        let nu = u0 * d / re;
        let tau = 0.5 + 3.0 * nu;
        let mut f = vec![0.0; nx * ny * nz * 19];
        for n in 0..nx * ny * nz {
            for i in 0..19 {
                f[n * 19 + i] = feq(i, 1.0, (u0, 0.0, 0.0));
            }
        }
        Lbm3 { nx, ny, nz, tau, u0, ftmp: f.clone(), f, body }
    }

    fn macros(&self, n: usize) -> (f64, (f64, f64, f64)) {
        let (mut rho, mut mx, mut my, mut mz) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..19 {
            let fi = self.f[n * 19 + i];
            rho += fi;
            mx += fi * E[i].0 as f64;
            my += fi * E[i].1 as f64;
            mz += fi * E[i].2 as f64;
        }
        if rho <= 0.0 { (1.0, (0.0, 0.0, 0.0)) } else { (rho, (mx / rho, my / rho, mz / rho)) }
    }

    pub fn step(&mut self) -> f64 {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let n = (z * ny + y) * nx + x;
                    if self.body.is_solid(x, y, z) { continue; }
                    let (rho, u) = self.macros(n);
                    for i in 0..19 {
                        let fi = self.f[n * 19 + i];
                        self.f[n * 19 + i] = fi - (fi - feq(i, rho, u)) / self.tau;
                    }
                }
            }
        }
        let mut drag = 0.0;
        for v in self.ftmp.iter_mut() { *v = 0.0; }
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let n = (z * ny + y) * nx + x;
                    if self.body.is_solid(x, y, z) { continue; }
                    for i in 0..19 {
                        let xn = x as i32 + E[i].0;
                        let yn = y as i32 + E[i].1;
                        let zn = z as i32 + E[i].2;
                        let fi = self.f[n * 19 + i];
                        if xn < 0 || xn >= nx as i32 || yn < 0 || yn >= ny as i32
                            || zn < 0 || zn >= nz as i32
                        {
                            self.ftmp[n * 19 + OPP[i]] += fi;
                            continue;
                        }
                        let (xu, yu, zu) = (xn as usize, yn as usize, zn as usize);
                        if self.body.is_solid(xu, yu, zu) {
                            self.ftmp[n * 19 + OPP[i]] += fi;
                            drag += 2.0 * E[i].0 as f64 * fi;
                        } else {
                            self.ftmp[((zu * ny + yu) * nx + xu) * 19 + i] += fi;
                        }
                    }
                }
            }
        }
        std::mem::swap(&mut self.f, &mut self.ftmp);
        self.inflow_outflow();
        drag
    }

    fn inflow_outflow(&mut self) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        for z in 0..nz {
            for y in 0..ny {
                let n = (z * ny + y) * nx; // x=0
                if !self.body.is_solid(0, y, z) {
                    for i in 0..19 { self.f[n * 19 + i] = feq(i, 1.0, (self.u0, 0.0, 0.0)); }
                }
                let no = (z * ny + y) * nx + (nx - 1);
                let np = (z * ny + y) * nx + (nx - 2);
                if !self.body.is_solid(nx - 1, y, z) {
                    for i in 0..19 { self.f[no * 19 + i] = self.f[np * 19 + i]; }
                }
            }
        }
    }

    /// Run `steps`; return the vehicle Cd = F_x / (½·ρ·u0²·A_frontal).
    pub fn run(&mut self, steps: usize) -> f64 {
        let tail = (steps / 5).max(1);
        let mut sum = 0.0;
        for s in 0..steps {
            let d = self.step();
            if s >= steps - tail { sum += d; }
        }
        let favg = sum / tail as f64;
        favg / (0.5 * self.u0 * self.u0 * self.body.frontal_cells as f64)
    }
}

/// Vehicle Cd of a 3D body at a Reynolds number.
pub fn vehicle_cd(body: Body3, re: f64, steps: usize) -> f64 {
    Lbm3::new(body, re).run(steps)
}
