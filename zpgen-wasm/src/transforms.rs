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

    // l1: source (origin) to zone plate point
    let l1 = norm2(r_o);

    // Image point at distance q0 from ZP center along optical axis
    let q_point = [
        p[0] + q0 * p_hat[0],
        p[1] + q0 * p_hat[1],
        p[2] + q0 * p_hat[2],
    ];

    // l2: zone plate point to image point (Euclidean)
    let l2_vec = [q_point[0] - r_o[0], q_point[1] - r_o[1], q_point[2] - r_o[2]];
    let l2 = norm2(&l2_vec);

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

/// Convert spatial frequency to zone plate XYZ coordinates (optimized version with pre-normalized normal)
///
/// This is an optimized version that assumes n_hat is already normalized,
/// avoiding redundant normalization in hot loops.
pub fn freq_to_zp_xyz_normalized(
    f: &[f64; 2],
    n_hat: &[f64; 3],
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

    // 2) Use pre-normalized plane normal (skip normalization)
    let n = n_hat;

    // 3) Choose lz sign so intersection is in front of the origin (t >= 0) if possible
    // Compute denom with current sign
    let mut denom = scalar_dot_product(&l, n);
    let numer = scalar_dot_product(p_on_plane, n);

    if denom == 0.0 {
        // Try flipping lz once; if still zero, it's parallel
        l[2] = -l[2];
        denom = scalar_dot_product(&l, n);
        if denom == 0.0 {
            return None;
        }
    }

    let mut t = numer / denom;

    // If you require forward intersection only (t >= 0), try flipping lz once to see if that helps
    if t < 0.0 {
        l[2] = -l[2];
        denom = scalar_dot_product(&l, n);
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

    #[test]
    fn test_xyz_to_opl_tilted_finite_conjugate_symmetry() {
        // Regression test: tilted ZP with finite conjugate should have
        // symmetric OPD on both sides of the tilt axis.
        // With the old plane-projection model, one side produced a spurious
        // linear OPD term causing grating artifacts.
        let lambda = 0.0135; // 13.5 nm in um
        let tilt = 6.0_f64.to_radians(); // 6 degree tilt about y-axis
        let dist = 100.0;

        // Tilted p vector (tilted in xz-plane)
        let p = [dist * tilt.sin(), 0.0, dist * tilt.cos()];
        // Finite conjugate: q0 = |p| (equal conjugates)
        let q0 = norm2(&p);

        // Tilted basis vectors
        let bx = [tilt.cos(), 0.0, -tilt.sin()];
        let by = [0.0, 1.0, 0.0];

        // Reference OPL at center
        let opl_center = xyz_to_opl(&p, &p, q0, lambda);

        // Test OPD at symmetric points along the tilt axis (bx direction)
        let r = 2.0; // 2 um offset
        let r_plus = [
            p[0] + r * bx[0],
            p[1] + r * bx[1],
            p[2] + r * bx[2],
        ];
        let r_minus = [
            p[0] - r * bx[0],
            p[1] - r * bx[1],
            p[2] - r * bx[2],
        ];

        let opd_plus = xyz_to_opl(&r_plus, &p, q0, lambda) - opl_center;
        let opd_minus = xyz_to_opl(&r_minus, &p, q0, lambda) - opl_center;

        // OPD should be approximately symmetric (both positive, similar magnitude)
        // The old bug caused one side to have a huge linear term
        let ratio = opd_plus / opd_minus;
        assert!(
            (ratio - 1.0).abs() < 0.1,
            "OPD should be approximately symmetric: opd_plus={}, opd_minus={}, ratio={}",
            opd_plus, opd_minus, ratio
        );

        // Both should be positive (longer path than on-axis)
        assert!(opd_plus > 0.0, "opd_plus should be positive: {}", opd_plus);
        assert!(opd_minus > 0.0, "opd_minus should be positive: {}", opd_minus);
    }

    #[test]
    fn test_xyz_to_opl_infinite_conjugate_backward_compat() {
        // With q0 >> |p|, the new point-based formula should give
        // results very close to the old plane-based formula.
        let lambda = 0.0135;
        let p = [0.0, 0.0, 100.0];
        let q0 = 100000.0; // quasi-infinite conjugate

        // Reference at center
        let opl_center = xyz_to_opl(&p, &p, q0, lambda);

        // Test point offset from center
        let r_o = [1.0, 0.0, 100.0];
        let opd = xyz_to_opl(&r_o, &p, q0, lambda) - opl_center;

        // For infinite conjugate, OPD ≈ r^2 / (2*f*lambda) where f ≈ p_norm
        // r = 1 um, f ≈ 100 um -> OPD ≈ 1 / (200 * 0.0135) ≈ 0.37 waves
        assert!(
            opd > 0.3 && opd < 0.5,
            "OPD for infinite conjugate should be ~0.37 waves, got {}",
            opd
        );
    }
}
