//! kami-cfd — clean-room D2Q9 lattice-Boltzmann fluid solver.
//!
//! The high-fidelity `:lbm` backend for the same CAE aero contract that
//! `aero-clj` answers with the reduced-order `:rom-buildup` build-up. Where the
//! build-up *correlates* drag from shape descriptors, this *resolves* the flow
//! field on a lattice and measures the drag on the body by momentum exchange
//! (Ladd). It returns a 2D **sectional** drag coefficient — i.e. relative shape
//! drag, useful for ranking afterbody/blockage choices; a full 3D vehicle Cd
//! needs the 3D extension (D3Q19/D3Q27), the same solver in another dimension.
//!
//! D2Q9, single-relaxation-time BGK. No external dependencies.
//!
//! The 3D extension (`d3`, D3Q19) yields a true vehicle Cd normalised by
//! frontal area — see [`d3::vehicle_cd`].

pub mod d3;
pub mod mesh;

/// D2Q9 lattice velocities.
const E: [(i32, i32); 9] = [
    (0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),
    (1, 1), (-1, 1), (-1, -1), (1, -1),
];
/// D2Q9 weights.
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
];
/// Opposite direction index (for bounce-back).
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// A 2D body profile placed in the channel, defined by a solid mask.
pub struct Body {
    pub nx: usize,
    pub ny: usize,
    pub solid: Vec<bool>,
    /// Body frontal height in lattice units (the reference length D).
    pub height: usize,
}

impl Body {
    /// A rectangular block of width `w`, height `h`, leading edge at `x0`,
    /// vertically centred — the bluff/squareback reference.
    pub fn block(nx: usize, ny: usize, x0: usize, w: usize, h: usize) -> Body {
        let mut solid = vec![false; nx * ny];
        let y0 = (ny - h) / 2;
        for x in x0..(x0 + w).min(nx) {
            for y in y0..(y0 + h).min(ny) {
                solid[y * nx + x] = true;
            }
        }
        Body { nx, ny, solid, height: h }
    }

    /// A tapered teardrop of the same frontal height `h`, length `len`, tail
    /// shrinking linearly to a point — the streamlined reference.
    pub fn teardrop(nx: usize, ny: usize, x0: usize, len: usize, h: usize) -> Body {
        let mut solid = vec![false; nx * ny];
        let cy = ny / 2;
        for k in 0..len {
            let x = x0 + k;
            if x >= nx { break; }
            // full height at the nose, tapering to ~0 at the tail
            let frac = 1.0 - (k as f64) / (len as f64);
            let half = ((h as f64) * 0.5 * frac).round() as usize;
            for y in cy.saturating_sub(half)..=(cy + half).min(ny - 1) {
                solid[y * nx + x] = true;
            }
        }
        Body { nx, ny, solid, height: h }
    }

    #[inline]
    fn is_solid(&self, x: usize, y: usize) -> bool {
        self.solid[y * self.nx + x]
    }
}

/// D2Q9 BGK lattice-Boltzmann solver over a channel with a body.
pub struct Lbm {
    nx: usize,
    ny: usize,
    tau: f64,
    u0: f64,
    f: Vec<f64>,
    ftmp: Vec<f64>,
    body: Body,
}

#[inline]
fn feq(i: usize, rho: f64, ux: f64, uy: f64) -> f64 {
    let (ex, ey) = E[i];
    let eu = ex as f64 * ux + ey as f64 * uy;
    let usq = ux * ux + uy * uy;
    W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq)
}

impl Lbm {
    /// Build a solver for a body at inflow speed set by `re` (Reynolds number
    /// on the body height). `u0` lattice velocity kept small for stability.
    pub fn new(body: Body, re: f64) -> Lbm {
        let (nx, ny) = (body.nx, body.ny);
        let u0 = 0.05;
        let nu = u0 * body.height as f64 / re; // ν = u0·D/Re
        let tau = 0.5 + 3.0 * nu;
        let mut f = vec![0.0; nx * ny * 9];
        for n in 0..nx * ny {
            for i in 0..9 {
                f[n * 9 + i] = feq(i, 1.0, u0, 0.0);
            }
        }
        Lbm { nx, ny, tau, u0, ftmp: f.clone(), f, body }
    }

    fn macros(&self, n: usize) -> (f64, f64, f64) {
        let mut rho = 0.0;
        let (mut mx, mut my) = (0.0, 0.0);
        for i in 0..9 {
            let fi = self.f[n * 9 + i];
            rho += fi;
            mx += fi * E[i].0 as f64;
            my += fi * E[i].1 as f64;
        }
        if rho <= 0.0 { (1.0, 0.0, 0.0) } else { (rho, mx / rho, my / rho) }
    }

    /// One collide + stream step; returns the x-drag on the body this step
    /// (momentum exchanged at bounce-back links).
    pub fn step(&mut self) -> f64 {
        let (nx, ny) = (self.nx, self.ny);
        // collision (BGK), in place
        for y in 0..ny {
            for x in 0..nx {
                let n = y * nx + x;
                if self.body.is_solid(x, y) { continue; }
                let (rho, ux, uy) = self.macros(n);
                for i in 0..9 {
                    let fi = self.f[n * 9 + i];
                    self.f[n * 9 + i] = fi - (fi - feq(i, rho, ux, uy)) / self.tau;
                }
            }
        }
        // streaming + bounce-back, accumulating drag
        let mut drag = 0.0;
        for v in self.ftmp.iter_mut() { *v = 0.0; }
        for y in 0..ny {
            for x in 0..nx {
                let n = y * nx + x;
                if self.body.is_solid(x, y) { continue; }
                for i in 0..9 {
                    let xn = x as i32 + E[i].0;
                    let yn = y as i32 + E[i].1;
                    let fi = self.f[n * 9 + i];
                    // walls (top/bottom) and out-of-range → bounce back to self
                    if yn < 0 || yn >= ny as i32 || xn < 0 || xn >= nx as i32 {
                        self.ftmp[n * 9 + OPP[i]] += fi;
                        continue;
                    }
                    let (xn, yn) = (xn as usize, yn as usize);
                    if self.body.is_solid(xn, yn) {
                        // bounce-back off the body; momentum exchange (Ladd)
                        self.ftmp[n * 9 + OPP[i]] += fi;
                        drag += 2.0 * E[i].0 as f64 * fi;
                    } else {
                        self.ftmp[(yn * nx + xn) * 9 + i] += fi;
                    }
                }
            }
        }
        std::mem::swap(&mut self.f, &mut self.ftmp);
        self.apply_boundaries();
        drag
    }

    /// Left inflow (equilibrium at u0), right outflow (zero-gradient).
    fn apply_boundaries(&mut self) {
        let (nx, ny) = (self.nx, self.ny);
        for y in 0..ny {
            let n = y * nx; // x = 0
            if !self.body.is_solid(0, y) {
                for i in 0..9 { self.f[n * 9 + i] = feq(i, 1.0, self.u0, 0.0); }
            }
            let no = y * nx + (nx - 1); // x = nx-1
            let np = y * nx + (nx - 2);
            if !self.body.is_solid(nx - 1, y) {
                for i in 0..9 { self.f[no * 9 + i] = self.f[np * 9 + i]; }
            }
        }
    }

    /// Run `steps` and return the sectional drag coefficient, averaged over the
    /// last 20% of steps. Cd = F_x / (½·ρ·u0²·D), ρ≈1.
    pub fn run(&mut self, steps: usize) -> f64 {
        let tail = (steps / 5).max(1);
        let mut sum = 0.0;
        for s in 0..steps {
            let d = self.step();
            if s >= steps - tail { sum += d; }
        }
        let favg = sum / tail as f64;
        favg / (0.5 * self.u0 * self.u0 * self.body.height as f64)
    }
}

/// Convenience: sectional Cd of a body at a Reynolds number.
pub fn sectional_cd(body: Body, re: f64, steps: usize) -> f64 {
    Lbm::new(body, re).run(steps)
}
