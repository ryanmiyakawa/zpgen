// Zernike polynomial functions
// Port from zpUtils.cpp

/// n Choose k function (binomial coefficient)
///
/// Port from zpUtils.cpp lines 11-18
fn n_choose_k(n: i32, k: i32) -> i32 {
    if k == 0 {
        return 1;
    }
    (n * n_choose_k(n - 1, k - 1)) / k
}

/// Retrieves the value of the Zernike Polynomial order J at the coordinate [r, th]
///
/// Uses Fringe/University of Arizona indexing convention
/// Port from zpUtils.cpp lines 21-77
pub fn zgenpt(j: i32, r: f64, th: f64) -> f64 {
    // Get dual index [n,m] from j:
    let mut smct = 2;
    let mut ct = 0;
    let mut num_in_smct = 3;
    let mut j_temp = j;

    while j_temp > 0 {
        smct += 2;
        j_temp -= num_in_smct;
        num_in_smct = smct + 1;
        ct += 1;
    }

    let mut n = 2 * ct;
    let mut m: i32 = 0;

    for q in 1..=j_temp.abs() {
        if q % 2 == 1 {
            n -= 1;
            m += 1;
        }
    }

    if j_temp.abs() % 2 == 1 {
        m *= -1;
    }

    // Compute value of radial function:
    let p = (n - m.abs()) / 2;
    let mut rv = 0.0;

    for k in 0..=p {
        let k_parity = if k % 2 == 1 { -1.0 } else { 1.0 };

        rv += (k_parity
            * n_choose_k(n - k, k) as f64
            * n_choose_k(n - 2 * k, p - k) as f64)
            * r.powi((n - 2 * k) as i32);
    }

    // Compute value of azimuthal function:
    let av = if m > 0 {
        (m as f64 * th).cos()
    } else if m < 0 {
        -(m as f64 * th).sin()
    } else {
        1.0
    };

    rv * av
}

/// Computes the phase of Zernike polynomial
///
/// Port from zpUtils.cpp lines 116-126
pub fn compute_zernike_phase(cx: f64, cy: f64, orders: &[f64], n_zerns: usize) -> f64 {
    let rn = (cx * cx + cy * cy).sqrt();
    let th = cy.atan2(cx);
    let mut zv = 0.0;

    for k in 0..n_zerns {
        zv += orders[n_zerns + k] * zgenpt(orders[k] as i32, rn, th);
    }

    zv
}

/// Returns phase term for the provided parameters
///
/// Port from zpUtils.cpp lines 128-137
pub fn get_phase_term(
    cx: f64,
    cy: f64,
    orders: &[f64],
    n_zerns: usize,
    zpc_phase: f64,
    zpcr1: f64,
    zpcr2: f64,
) -> f64 {
    let r = (cx * cx + cy * cy).sqrt();
    let mut ph = 0.0;

    if r <= zpcr1 && r >= zpcr2 {
        ph += zpc_phase / (2.0 * std::f64::consts::PI);
    }

    ph + compute_zernike_phase(cx, cy, orders, n_zerns)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_choose_k() {
        assert_eq!(n_choose_k(5, 0), 1);
        assert_eq!(n_choose_k(5, 1), 5);
        assert_eq!(n_choose_k(5, 2), 10);
        assert_eq!(n_choose_k(5, 3), 10);
    }

    #[test]
    fn test_zgenpt_piston() {
        // Zernike index 0 should be piston (constant = 1)
        let val = zgenpt(0, 0.5, 0.0);
        assert!((val - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zgenpt_tilt() {
        // Test that Zernike polynomials are computed
        let val = zgenpt(1, 0.5, std::f64::consts::PI / 4.0);
        assert!(val.is_finite());
    }
}
