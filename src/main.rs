//! kami-cfd CLI — compute the sectional Cd of a built-in body profile.
//!
//! Usage: kami-cfd [block|teardrop] [re] [steps]
//! Emits an EDN-ish line the `:lbm` cae-solver method can read (subprocess /
//! FFI wiring is the next step; the solver contract slot already exists).

use kami_cfd::d3::{vehicle_cd, Body3};
use kami_cfd::{sectional_cd, Body};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let shape = args.get(1).map(String::as_str).unwrap_or("block");
    let re: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(100.0);
    let steps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3000);

    // 3D mode: true vehicle Cd (frontal-area-normalised). 2D mode: sectional.
    if shape == "box3d" || shape == "fastback3d" {
        let (nx, ny, nz) = (96, 40, 28);
        let body = if shape == "fastback3d" {
            Body3::fastback_car(nx, ny, nz, 16, 44, 20, 16)
        } else {
            Body3::box_car(nx, ny, nz, 16, 44, 20, 16)
        };
        let cd = vehicle_cd(body, re, steps.max(1500));
        println!(
            "{{:solver :lbm :dim 3 :shape :{} :re {} :steps {} :vehicle-cd {:.4}}}",
            shape, re, steps.max(1500), cd
        );
        return;
    }

    let (nx, ny) = (240, 80);
    let body = match shape {
        "teardrop" => Body::teardrop(nx, ny, 40, 60, 24),
        _ => Body::block(nx, ny, 40, 24, 24),
    };
    let cd = sectional_cd(body, re, steps);
    println!(
        "{{:solver :lbm :dim 2 :shape :{} :re {} :steps {} :sectional-cd {:.4}}}",
        shape, re, steps, cd
    );
}
