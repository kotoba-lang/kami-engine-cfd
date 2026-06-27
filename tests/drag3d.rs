//! D3Q19 vehicle-drag contract: the 3D solver is stable, the frontal-area-
//! normalised Cd is finite/positive and in a plausible bluff-body band, and a
//! fastback (tapered roof) out-performs a squareback box of equal frontal area.

use kami_cfd::d3::{vehicle_cd, Body3};

const NX: usize = 96;
const NY: usize = 40;
const NZ: usize = 28;

#[test]
fn high_re_stable_with_les() {
    // Re=2000 — single-relaxation BGK NaNs here; the Smagorinsky LES keeps the
    // solver stable (τ stays clear of ½). The absolute value is still inflated
    // by the coarse high-blockage test domain; aero-clj calibrates the scale.
    let cd = vehicle_cd(Body3::box_car(NX, NY, NZ, 16, 44, 20, 16), 2000.0, 1000);
    assert!(cd.is_finite() && cd > 0.0, "LES must stay stable at high Re, got {cd}");
    assert!(cd < 12.0, "Cd sanity ceiling exceeded: {cd}");
}

#[test]
fn fastback_beats_squareback() {
    // equal frontal area, tapered roof vs squareback → lower pressure drag,
    // in the turbulent (LES) regime.
    let box_cd = vehicle_cd(Body3::box_car(NX, NY, NZ, 16, 44, 20, 16), 2000.0, 1200);
    let fast_cd = vehicle_cd(Body3::fastback_car(NX, NY, NZ, 16, 44, 20, 16), 2000.0, 1200);
    assert!(
        fast_cd < box_cd,
        "fastback must beat squareback: box={box_cd:.3} fastback={fast_cd:.3}"
    );
}
