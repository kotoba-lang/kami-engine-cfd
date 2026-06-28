//! Mesh import + parametric geometry — real vehicle surfaces voxelised into the
//! D3Q19 solver, replacing the box/fastback primitives. Binary STL parser plus
//! an Ahmed-body generator (the canonical automotive CFD validation bluff body:
//! a block with a slanted rear). Winding is irrelevant downstream — the
//! voxeliser uses ray-parity, which only counts surface crossings.

/// A triangle as three (x,y,z) vertices.
pub type Tri = [[f64; 3]; 3];

/// Parse a binary STL (80-byte header, u32 triangle count, 50 bytes/triangle:
/// 3 floats normal + 3×3 floats verts + 2 attr bytes).
pub fn parse_stl(bytes: &[u8]) -> Vec<Tri> {
    if bytes.len() < 84 {
        return vec![];
    }
    let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    let rdf = |o: usize| f32::from_le_bytes([bytes[o], bytes[o + 1], bytes[o + 2], bytes[o + 3]]) as f64;
    let mut tris = Vec::with_capacity(n);
    let mut off = 84;
    for _ in 0..n {
        if off + 50 > bytes.len() {
            break;
        }
        let mut t = [[0.0; 3]; 3];
        for v in 0..3 {
            let b = off + 12 + v * 12; // skip the 12-byte normal
            t[v] = [rdf(b), rdf(b + 4), rdf(b + 8)];
        }
        tris.push(t);
        off += 50;
    }
    tris
}

/// Extrude a side profile (CCW list of (x,z) points) along y from 0..width into
/// a closed triangle soup.
fn prism(profile: &[(f64, f64)], width: f64) -> Vec<Tri> {
    let n = profile.len();
    let p = |i: usize, y: f64| [profile[i].0, y, profile[i].1];
    let mut tris = vec![];
    for i in 1..n - 1 {
        tris.push([p(0, 0.0), p(i, 0.0), p(i + 1, 0.0)]);     // near cap
        tris.push([p(0, width), p(i + 1, width), p(i, width)]); // far cap
    }
    for i in 0..n {
        let j = (i + 1) % n;
        tris.push([p(i, 0.0), p(i, width), p(j, 0.0)]);
        tris.push([p(j, 0.0), p(i, width), p(j, width)]);
    }
    tris
}

/// Ahmed-like body: length 10, height 3, width 4, with a rear-roof slant.
/// `slant_deg >= 90` → squareback (vertical rear). Sits on z=0.
pub fn ahmed_body(slant_deg: f64) -> Vec<Tri> {
    let (l, h, w) = (10.0, 3.0, 4.0);
    let profile: Vec<(f64, f64)> = if slant_deg >= 90.0 {
        vec![(0.0, 0.0), (l, 0.0), (l, h), (0.0, h)]
    } else {
        let slant_x = 3.0;
        let hr = (h - slant_x * slant_deg.to_radians().tan()).max(0.2);
        vec![(0.0, 0.0), (l, 0.0), (l, hr), (l - slant_x, h), (0.0, h)]
    };
    prism(&profile, w)
}
