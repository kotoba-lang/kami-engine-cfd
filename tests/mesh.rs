//! Real-geometry pipeline: a triangle mesh voxelises into a solid body and runs
//! through the D3Q19 LES solver. Validation geometry = the Ahmed body (canonical
//! automotive bluff body); a rear slant must reduce drag vs a squareback.

use kami_cfd::d3::{vehicle_cd, Body3, Lbm3};
use kami_cfd::mesh;

const NX: usize = 140;
const NY: usize = 64;
const NZ: usize = 40;

#[test]
fn mesh_voxelises_to_a_solid() {
    let tris = mesh::ahmed_body(90.0); // squareback box
    let body = Body3::from_triangles(NX, NY, NZ, &tris, 48.0, 20.0);
    let filled = body.solid.iter().filter(|&&s| s).count();
    assert!(filled > 1000, "voxelised body should be substantial, got {filled} cells");
    assert!(body.frontal_cells > 100, "frontal area too small: {}", body.frontal_cells);
}

#[test]
fn ahmed_slant_reduces_drag() {
    // a 25° rear slant must out-perform a vertical squareback (same box).
    let square = Body3::from_triangles(NX, NY, NZ, &mesh::ahmed_body(90.0), 48.0, 20.0);
    let slant = Body3::from_triangles(NX, NY, NZ, &mesh::ahmed_body(25.0), 48.0, 20.0);
    let sq_cd = Lbm3::new(square, 2000.0).run(1200);
    let sl_cd = Lbm3::new(slant, 2000.0).run(1200);
    assert!(sq_cd.is_finite() && sl_cd.is_finite());
    assert!(
        sl_cd < sq_cd,
        "25° slant must beat squareback: square={sq_cd:.3} slant={sl_cd:.3}"
    );
}

#[test]
fn stl_roundtrips() {
    // a tiny hand-built binary STL (1 triangle) parses back to one triangle.
    let mut b = vec![0u8; 84 + 50];
    b[80] = 1; // count = 1 (little-endian u32)
    let put = |buf: &mut [u8], o: usize, f: f32| buf[o..o + 4].copy_from_slice(&f.to_le_bytes());
    // verts at offset 84+12 .. (skip normal)
    let mut o = 84 + 12;
    for &c in &[0.0f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0] {
        put(&mut b, o, c);
        o += 4;
    }
    let tris = mesh::parse_stl(&b);
    assert_eq!(tris.len(), 1);
    assert_eq!(tris[0][1], [1.0, 0.0, 0.0]);
    let _ = vehicle_cd; // silence unused import on some configs
}
