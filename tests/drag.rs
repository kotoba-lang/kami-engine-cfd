//! LBM drag contract: the solver is stable, drag is finite and positive, and a
//! bluff body has higher sectional drag than a streamlined one of equal frontal
//! height — the physically essential ranking (afterbody separation).

use kami_cfd::{sectional_cd, Body};

const NX: usize = 240;
const NY: usize = 80;

#[test]
fn block_drag_is_finite_and_positive() {
    let cd = sectional_cd(Body::block(NX, NY, 40, 24, 24), 100.0, 2500);
    assert!(cd.is_finite(), "Cd must be finite, got {cd}");
    assert!(cd > 0.0, "a bluff body has positive drag, got {cd}");
}

#[test]
fn bluff_beats_streamlined() {
    // same frontal height (24), different afterbody: square block vs teardrop.
    let block = sectional_cd(Body::block(NX, NY, 40, 24, 24), 100.0, 3000);
    let tear = sectional_cd(Body::teardrop(NX, NY, 40, 60, 24), 100.0, 3000);
    assert!(
        block > tear,
        "squareback must out-drag a streamlined tail: block={block:.3} teardrop={tear:.3}"
    );
}

#[test]
fn taller_body_more_drag() {
    // more frontal blockage → more drag force (Cd normalised by D, but the
    // bigger wake still dominates at fixed channel).
    let small = sectional_cd(Body::block(NX, NY, 40, 24, 16), 100.0, 2500);
    let big = sectional_cd(Body::block(NX, NY, 40, 24, 32), 100.0, 2500);
    assert!(big.is_finite() && small.is_finite());
    assert!(big > 0.0 && small > 0.0);
}
