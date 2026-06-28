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
/// Specular reflection across a y-wall: (ex,ey,ez) → (ex,-ey,ez). Free-slip.
const REFLECT_Y: [usize; 19] =
    [0, 1, 2, 4, 3, 5, 6, 9, 10, 7, 8, 11, 12, 13, 14, 18, 17, 16, 15];
/// Specular reflection across a z-wall: (ex,ey,ez) → (ex,ey,-ez). Free-slip.
const REFLECT_Z: [usize; 19] =
    [0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 13, 14, 11, 12, 17, 18, 15, 16];

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

    /// Voxelise a triangle mesh (REAL geometry, e.g. an STL or ahmed_body) into
    /// the domain. The mesh is uniformly scaled so its length maps to
    /// `len_cells`, centred in y, sitting on the ground (z=1). Solid filling is
    /// per-(y,z)-column ray parity: cast a +x ray through each column, sort the
    /// surface crossings, fill between entry/exit pairs (watertight meshes).
    pub fn from_triangles(
        nx: usize, ny: usize, nz: usize,
        tris: &[crate::mesh::Tri], len_cells: f64, x0: f64,
    ) -> Body3 {
        // mesh bounding box
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for t in tris {
            for v in t {
                for a in 0..3 {
                    lo[a] = lo[a].min(v[a]);
                    hi[a] = hi[a].max(v[a]);
                }
            }
        }
        let scale = len_cells / (hi[0] - lo[0]).max(1e-9);
        let ymid = 0.5 * (lo[1] + hi[1]);
        let tf = |v: [f64; 3]| -> [f64; 3] {
            [x0 + (v[0] - lo[0]) * scale,
             ny as f64 / 2.0 + (v[1] - ymid) * scale,
             1.0 + (v[2] - lo[2]) * scale]
        };
        let lt: Vec<crate::mesh::Tri> = tris.iter().map(|t| [tf(t[0]), tf(t[1]), tf(t[2])]).collect();
        let mut solid = vec![false; nx * ny * nz];
        for z in 0..nz {
            for y in 0..ny {
                let (uc, vc) = (y as f64 + 0.5, z as f64 + 0.5); // column centre (y,z)
                let mut xs: Vec<f64> = vec![];
                for t in &lt {
                    if let Some(x) = ray_x_cross(t, uc, vc) {
                        xs.push(x);
                    }
                }
                if xs.len() < 2 { continue; }
                xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mut k = 0;
                while k + 1 < xs.len() {
                    let a = xs[k].ceil().max(0.0) as i64;
                    let b = xs[k + 1].floor().min((nx - 1) as f64) as i64;
                    for x in a..=b {
                        if x >= 0 && (x as usize) < nx {
                            solid[(z * ny + y) * nx + x as usize] = true;
                        }
                    }
                    k += 2;
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

/// Intersection x of a +x ray at (y=uc, z=vc) with a triangle, via barycentric
/// containment in the y-z projection. None if the ray misses or the triangle is
/// edge-on to the projection.
fn ray_x_cross(t: &crate::mesh::Tri, uc: f64, vc: f64) -> Option<f64> {
    let (u0, v0) = (t[0][1], t[0][2]);
    let (u1, v1) = (t[1][1], t[1][2]);
    let (u2, v2) = (t[2][1], t[2][2]);
    let det = (v1 - v2) * (u0 - u2) + (u2 - u1) * (v0 - v2);
    if det.abs() < 1e-9 {
        return None;
    }
    let a = ((v1 - v2) * (uc - u2) + (u2 - u1) * (vc - v2)) / det;
    let b = ((v2 - v0) * (uc - u2) + (u0 - u2) * (vc - v2)) / det;
    let c = 1.0 - a - b;
    if a < -1e-9 || b < -1e-9 || c < -1e-9 {
        return None;
    }
    Some(a * t[0][0] + b * t[1][0] + c * t[2][0])
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

    /// BGK collision with a Smagorinsky LES subgrid model — per-cell
    /// independent, split across threads (std::thread::scope, zero-dep). The
    /// local strain rate is read from the non-equilibrium momentum-flux tensor
    /// Q_ab = Σ_i e_ia·e_ib·(f_i−feq_i); a turbulent eddy viscosity
    /// ν_t = (Cs·Δ)²·S̄ raises τ where the flow is sheared. This keeps τ clear
    /// of ½ at high Re (stability) AND models subgrid turbulence — the physical
    /// requirement for the drag to fall toward the real attached-flow value,
    /// which single-relaxation BGK alone cannot do.
    fn collide_parallel(&mut self) {
        let ncells = self.nx * self.ny * self.nz;
        let solid = &self.body.solid;
        let tau0 = self.tau;
        let cs = 0.16; // Smagorinsky constant
        let nthreads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        let per = (ncells + nthreads - 1) / nthreads;
        std::thread::scope(|s| {
            for (t, fc) in self.f.chunks_mut(per * 19).enumerate() {
                let base = t * per;
                s.spawn(move || {
                    let cells = fc.len() / 19;
                    for c in 0..cells {
                        let n = base + c;
                        if solid[n] { continue; }
                        let (mut rho, mut mx, mut my, mut mz) = (0.0, 0.0, 0.0, 0.0);
                        for i in 0..19 {
                            let fi = fc[c * 19 + i];
                            rho += fi;
                            mx += fi * E[i].0 as f64;
                            my += fi * E[i].1 as f64;
                            mz += fi * E[i].2 as f64;
                        }
                        let u = if rho <= 0.0 { (0.0, 0.0, 0.0) } else { (mx / rho, my / rho, mz / rho) };
                        let rho = if rho <= 0.0 { 1.0 } else { rho };
                        // non-equilibrium momentum-flux tensor Q_ab
                        let (mut qxx, mut qyy, mut qzz, mut qxy, mut qxz, mut qyz) =
                            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                        let mut feqs = [0.0f64; 19];
                        for i in 0..19 {
                            feqs[i] = feq(i, rho, u);
                            let neq = fc[c * 19 + i] - feqs[i];
                            let (ex, ey, ez) = (E[i].0 as f64, E[i].1 as f64, E[i].2 as f64);
                            qxx += ex * ex * neq; qyy += ey * ey * neq; qzz += ez * ez * neq;
                            qxy += ex * ey * neq; qxz += ex * ez * neq; qyz += ey * ez * neq;
                        }
                        let qnorm = (2.0 * (qxx * qxx + qyy * qyy + qzz * qzz)
                            + 4.0 * (qxy * qxy + qxz * qxz + qyz * qyz)).sqrt();
                        // total relaxation time incl. Smagorinsky eddy viscosity
                        let tau_t = 0.5 * ((tau0 * tau0 + 18.0 * 1.4142135623730951 * cs * cs * qnorm / rho).sqrt() - tau0);
                        let omega = 1.0 / (tau0 + tau_t);
                        for i in 0..19 {
                            fc[c * 19 + i] -= omega * (fc[c * 19 + i] - feqs[i]);
                        }
                    }
                });
            }
        });
    }

    /// Streaming + drag, PULL form (gather) so it parallelises with no write
    /// conflicts: each cell writes only its own ftmp, reading upstream f. Wall/
    /// solid upstream → halfway bounce-back. Drag (momentum exchange) is a
    /// read-only reduction over fluid-solid links. Returns total x-drag.
    fn stream_drag_parallel(&mut self) -> f64 {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let solid = &self.body.solid;
        let f = &self.f;
        let ncells = nx * ny * nz;
        let nthreads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        let per = (ncells + nthreads - 1) / nthreads;
        let drags: Vec<f64> = std::thread::scope(|s| {
            let handles: Vec<_> = self.ftmp.chunks_mut(per * 19).enumerate().map(|(t, fc)| {
                let base = t * per;
                s.spawn(move || {
                    let cells = fc.len() / 19;
                    let mut drag = 0.0;
                    for c in 0..cells {
                        let n = base + c;
                        let z = n / (nx * ny);
                        let y = (n / nx) % ny;
                        let x = n % nx;
                        if solid[n] { for j in 0..19 { fc[c * 19 + j] = 0.0; } continue; }
                        for j in 0..19 {
                            let sx = x as i32 - E[j].0;
                            let sy = y as i32 - E[j].1;
                            let sz = z as i32 - E[j].2;
                            let val = if sz < 0 {
                                f[n * 19 + OPP[j]]                       // ground: no-slip road
                            } else if sy < 0 || sy >= ny as i32 {
                                f[n * 19 + REFLECT_Y[j]]                 // sides: free-slip far-field
                            } else if sz >= nz as i32 {
                                f[n * 19 + REFLECT_Z[j]]                 // top: free-slip far-field
                            } else if sx < 0 || sx >= nx as i32 {
                                f[n * 19 + OPP[j]]                       // x ends: set by inflow/outflow
                            } else {
                                let sn = ((sz as usize * ny) + sy as usize) * nx + sx as usize;
                                if solid[sn] { f[n * 19 + OPP[j]] }      // body bounce-back
                                else { f[sn * 19 + j] }                  // pull from upstream
                            };
                            fc[c * 19 + j] = val;
                        }
                        // momentum exchange on links pointing into the body
                        for i in 0..19 {
                            let dx = x as i32 + E[i].0;
                            let dy = y as i32 + E[i].1;
                            let dz = z as i32 + E[i].2;
                            if dx >= 0 && dx < nx as i32 && dy >= 0 && dy < ny as i32
                                && dz >= 0 && dz < nz as i32
                            {
                                let dn = ((dz as usize * ny) + dy as usize) * nx + dx as usize;
                                if solid[dn] { drag += 2.0 * E[i].0 as f64 * f[n * 19 + i]; }
                            }
                        }
                    }
                    drag
                })
            }).collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        std::mem::swap(&mut self.f, &mut self.ftmp);
        self.inflow_outflow();
        drags.iter().sum()
    }

    pub fn step(&mut self) -> f64 {
        self.collide_parallel();
        self.stream_drag_parallel()
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
