//! kami-cfd CLI — compute the sectional Cd of a built-in body profile.
//!
//! Usage: kami-cfd [block|teardrop] [re] [steps]
//! Emits an EDN-ish line the `:lbm` cae-solver method can read (subprocess /
//! FFI wiring is the next step; the solver contract slot already exists).

use kami_cfd::{sectional_cd, Body};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let shape = args.get(1).map(String::as_str).unwrap_or("block");
    let re: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(100.0);
    let steps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3000);

    let (nx, ny) = (240, 80);
    let body = match shape {
        "teardrop" => Body::teardrop(nx, ny, 40, 60, 24),
        _ => Body::block(nx, ny, 40, 24, 24),
    };
    let cd = sectional_cd(body, re, steps);
    println!(
        "{{:solver :lbm :shape :{} :re {} :steps {} :sectional-cd {:.4}}}",
        shape, re, steps, cd
    );
}
