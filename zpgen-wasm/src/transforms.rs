// Coordinate transformation utilities
// Port from zpUtils.cpp

/// Compute 2-norm of a 3D vector
#[inline]
pub fn norm2(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Compute scalar dot product of two 3D vectors
#[inline]
pub fn scalar_dot_product(v1: &[f64; 3], v2: &[f64; 3]) -> f64 {
    v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
}

/// Normalize a 3D vector in-place
pub fn norm_vector(v: &mut [f64; 3]) {
    let norm = norm2(v);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
}

/// Compute cross product of two 3D vectors
///
/// Port from zpUtils.cpp lines 505-509
pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Convert zone plate (Ux, Uy) coordinates to XYZ coordinates
///
/// Port from zpUtils.cpp lines 627-634
pub fn zpuxuy_to_xyz(u: &[f64; 3], p0: &[f64; 3], bx: &[f64; 3], by: &[f64; 3]) -> [f64; 3] {
    // Normalization of basis vectors is skipped because assuming they are normalized

    // Convert U to Cartesian and write in terms of p0 origin
    [
        p0[0] + u[0] * bx[0] + u[1] * by[0],
        p0[1] + u[0] * bx[1] + u[1] * by[1],
        p0[2] + u[0] * bx[2] + u[1] * by[2],
    ]
}

/// Convert XYZ coordinates to zone plate (Ux, Uy, Uz) coordinates
///
/// Port from zpUtils.cpp lines 611-625
pub fn zpxyz_to_uxuy(r: &[f64; 3], p: &[f64; 3], bx: &[f64; 3], by: &[f64; 3]) -> [f64; 3] {
    // Normalization of basis vectors is skipped because assuming they are normalized

    // Point in Cartesian coordinates wrt p0
    let rp = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];

    // compute bz from bx and by
    let bz = cross_product(bx, by);

    // Project onto basis
    [
        scalar_dot_product(&rp, bx),
        scalar_dot_product(&rp, by),
        scalar_dot_product(&rp, &bz),
    ]
}

/// Compute optical path length in wavelengths
///
/// Computes the path length from origin to r_o, then to image point at distance q0
/// Handles both virtual and real sources and objects
///
/// Port from zpUtils.cpp lines 645-690
pub fn xyz_to_opl(r_o: &[f64; 3], p: &[f64; 3], q0: f64, lambda: f64) -> f64 {
    let p_norm = norm2(p);
    let p_hat = [p[0] / p_norm, p[1] / p_norm, p[2] / p_norm];

    // distance from object (origin) to point r_o
    let l1 = norm2(r_o);

    // projection of r_o onto the optical axis
    let proj = r_o[0] * p_hat[0] + r_o[1] * p_hat[1] + r_o[2] * p_hat[2];

    // how much farther along p we need to go to reach the image plane
    let delta = q0 - proj;

    // compute image point r_i = r_o + delta * p_hat
    let r_i = [
        r_o[0] + delta * p_hat[0],
        r_o[1] + delta * p_hat[1],
        r_o[2] + delta * p_hat[2],
    ];

    // distance from r_o to r_i
    let l2_vec = [r_i[0] - r_o[0], r_i[1] - r_o[1], r_i[2] - r_o[2]];
    let l2 = norm2(&l2_vec);

    // return total optical path length in waves
    (l1 + l2) / lambda
}

/// Convert ZP coordinates to k-vector (spatial frequency)
///
/// Port from zpUtils.cpp lines 605-609
pub fn zpcoord_to_kvector(r: &[f64; 3]) -> [f64; 2] {
    let r_norm = norm2(r);
    [r[0] / r_norm, r[1] / r_norm]
}

/// Convert spatial frequency to zone plate XYZ coordinates
///
/// Port from zpUtils.cpp lines 522-603
pub fn freq_to_zp_xyz(
    f: &[f64; 2],
    n_in: &[f64; 3],
    p_on_plane: &[f64; 3],
    lambda: f64,
) -> Option<[f64; 3]> {
    let fx = f[0];
    let fy = f[1];

    // 1) Build unit direction l from spatial frequencies
    let s2 = (lambda * fx) * (lambda * fx) + (lambda * fy) * (lambda * fy);
    if s2 > 1.0 {
        // Evanescent: no real propagation direction
        return None;
    }

    let lz = (1.0 - s2).max(0.0).sqrt();
    let mut l = [lambda * fx, lambda * fy, lz]; // already unit length

    // 2) Normalize plane normal (without mutating input)
    let mut n = *n_in;
    let nmag = norm2(&n);
    if nmag == 0.0 {
        return None;
    }
    n[0] /= nmag;
    n[1] /= nmag;
    n[2] /= nmag;

    // 3) Choose lz sign so intersection is in front of the origin (t >= 0) if possible
    // Compute denom with current sign
    let mut denom = scalar_dot_product(&l, &n);
    let numer = scalar_dot_product(p_on_plane, &n);

    if denom == 0.0 {
        // Try flipping lz once; if still zero, it's parallel
        l[2] = -l[2];
        denom = scalar_dot_product(&l, &n);
        if denom == 0.0 {
            return None;
        }
    }

    let mut t = numer / denom;

    // If you require forward intersection only (t >= 0), try flipping lz once to see if that helps
    if t < 0.0 {
        l[2] = -l[2];
        denom = scalar_dot_product(&l, &n);
        if denom == 0.0 {
            return None;
        }
        t = numer / denom;
        if t < 0.0 {
            // Plane is behind the chosen propagation direction
            return None;
        }
    }

    // 4) Intersection point
    Some([t * l[0], t * l[1], t * l[2]])
}

/// Composite function: ZP (R, theta) to pupil coordinates (cx, cy)
///
/// Port from zpUtils.cpp lines 693-717
pub fn zp_rth_to_pc_cxcy(
    r: f64,
    th: f64,
    k_0: &[f64; 2],
    p: &[f64; 3],
    bx: &[f64; 3],
    by: &[f64; 3],
    lambda: f64,
    na: f64,
    anamorphic_azimuth: f64,
) -> [f64; 2] {
    let ux = r * th.cos();
    let uy = r * th.sin();
    let u = [ux, uy, 0.0];

    let r_xyz = zpuxuy_to_xyz(&u, p, bx, by);
    let k = zpcoord_to_kvector(&r_xyz);

    let cxy = [(k[0] - k_0[0]) / na, (k[1] - k_0[1]) / na];

    // Determine in-plane angle of k by taking atan2(k[1],k[0]):
    let k_azi = anamorphic_azimuth; // -atan2(k_0[1], k_0[0]);

    // Rotate C by k_azi:
    [
        k_azi.cos() * cxy[0] - k_azi.sin() * cxy[1],
        k_azi.sin() * cxy[0] + k_azi.cos() * cxy[1],
    ]
}

/// Convert frequency to zone plate Ux, Uy coordinates
///
/// Port from zpUtils.cpp lines 720-728
pub fn freq_to_zp_uxuy(
    fq: &[f64; 2],
    n: &[f64; 3],
    p: &[f64; 3],
    bx: &[f64; 3],
    by: &[f64; 3],
    lambda: f64,
) -> Option<[f64; 3]> {
    // F -> r
    let r = freq_to_zp_xyz(fq, n, p, lambda)?;

    // r -> U
    Some(zpxyz_to_uxuy(&r, p, bx, by))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_norm2() {
        let v = [3.0, 4.0, 0.0];
        assert_eq!(norm2(&v), 5.0);
    }

    #[test]
    fn test_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross_product(&a, &b);
        assert_eq!(c, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn test_zpuxuy_to_xyz_identity() {
        // Test with identity basis vectors and origin
        let u = [1.0, 2.0, 0.0];
        let p = [0.0, 0.0, 0.0];
        let bx = [1.0, 0.0, 0.0];
        let by = [0.0, 1.0, 0.0];

        let xyz = zpuxuy_to_xyz(&u, &p, &bx, &by);
        assert_eq!(xyz, [1.0, 2.0, 0.0]);
    }
}
