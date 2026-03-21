// Numerical solver for zone plate radii
// Port of objectiveFn and secantSolve from zpGenHolo.cpp

use crate::transforms::{zpuxuy_to_xyz, xyz_to_opl};

/// Objective function for zone plate radius solving
/// Returns the optical path difference for a given radius
///
/// Port of objectiveFn() from zpGenHolo.cpp lines 26-46
pub fn objective_fn(
    r: f64,
    th: f64,
    n: f64,
    p: &[f64; 3],
    q0: f64,
    bx: &[f64; 3],
    by: &[f64; 3],
    phase: f64,
    lambda: f64,
) -> f64 {
    // Compute xyz space coordinates from (r,th)
    // U are the bx,by coordinates (in-plane coordinates)
    let u = [r * th.cos(), r * th.sin(), 0.0];

    // Convert these back to xyz using basis vectors
    let r_xyz = zpuxuy_to_xyz(&u, p, bx, by);

    // Compute OPL in waves
    let opd = xyz_to_opl(&r_xyz, p, q0, lambda) - xyz_to_opl(p, p, q0, lambda);
    let zp_term = -n / 2.0;

    opd + zp_term + phase
}

/// Secant method solver for zone plate radii
///
/// Port of secantSolve() from zpGenHolo.cpp lines 48-74
pub fn secant_solve(
    dr_guess: f64,
    th: f64,
    n: f64,
    p: &[f64; 3],
    q0: f64,
    bx: &[f64; 3],
    by: &[f64; 3],
    phase: f64,
    lambda: f64,
) -> Result<f64, String> {
    const TOL_X: f64 = 0.00001;
    const MAX_ITER: usize = 50;

    let mut r1 = dr_guess;
    let mut r2 = r1 * 1.02;

    // Debug log initial guess
    // #[cfg(target_arch = "wasm32")]
    // if phase.abs() > 0.8 || dr_guess < 0.0 {
    //     web_sys::console::log_1(&format!(
    //         "Secant starting: theta={:.4}, phase={:.4}, initial_guess={:.6}, n={}",
    //         th, phase, dr_guess, n
    //     ).into());
    // }

    for iter in 0..MAX_ITER {
        let fr1 = objective_fn(r1, th, n, p, q0, bx, by, phase, lambda);
        let fr2 = objective_fn(r2, th, n, p, q0, bx, by, phase, lambda);

        let denominator = fr1 - fr2;

        // Check for problematic denominator
        if denominator.abs() < 1e-15 {
            return Err(format!(
                "SECANT DENOMINATOR TOO SMALL at theta={}, iter={}, fr1={}, fr2={}, r1={}, r2={}",
                th, iter, fr1, fr2, r1, r2
            ));
        }

        let r0 = r1 - fr1 * (r1 - r2) / denominator;

        // Check for invalid radius
        if !r0.is_finite() {
            return Err(format!(
                "SECANT PRODUCED NON-FINITE RADIUS at theta={}, iter={}, r0={}, fr1={}, fr2={}",
                th, iter, r0, fr1, fr2
            ));
        }

        if r0 < 0.0 {
            return Err(format!(
                "SECANT PRODUCED NEGATIVE RADIUS at theta={}, iter={}, r0={}, fr1={}, fr2={}",
                th, iter, r0, fr1, fr2
            ));
        }

        // Debug logging for spiral phase at problematic angles
        // #[cfg(target_arch = "wasm32")]
        // if phase > 0.4 && phase < 0.6 && iter < 5 {
        //     web_sys::console::log_1(&format!(
        //         "Secant iter {}: theta={:.4}, phase={:.4}, r1={:.6}, fr1={:.6}, r2={:.6}, fr2={:.6}, r0={:.6}",
        //         iter, th, phase, r1, fr1, r2, fr2, r0
        //     ).into());
        // }

        // Check convergence
        if (r0 - r1).abs() < TOL_X {
            return Ok(r0);
        }

        // Set new guesses
        r2 = r1;
        r1 = r0;
    }

    Err(format!(
        "MAXIMUM ITERATIONS REACHED at theta={}, phase={}, final_r={}, final_f={}",
        th, phase, r1, objective_fn(r1, th, n, p, q0, bx, by, phase, lambda)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secant_solve_simple() {
        // Test simple zone plate: p=[0,0,100], q=100000, NA=0.02
        let p = [0.0, 0.0, 100.0];
        let q = 100000.0;
        let bx = [1.0, 0.0, 0.0];
        let by = [0.0, 1.0, 0.0];
        let lambda = 0.0135; // 13.5 nm in um

        // Zone 10: r ≈ sqrt(n * lambda * f) where f = p*q/(p+q) ≈ 99.9
        // r ≈ sqrt(10 * 0.0135 * 99.9) ≈ 3.67 um
        let result = secant_solve(4.0, 0.0, 10.0, &p, q, &bx, &by, 0.0, lambda);

        assert!(result.is_ok());
        let radius = result.unwrap();
        assert!(radius > 2.0 && radius < 6.0, "Radius {} should be between 2 and 6 um", radius);
    }

    #[test]
    fn test_secant_solve_tilted_finite_conjugate() {
        // Regression test: tilted ZP with finite conjugate should produce
        // comparable radii at theta=0 and theta=PI (opposite sides of tilt).
        // The old bug caused one side to behave like a grating.
        use crate::transforms::norm2;

        let tilt = 6.0_f64.to_radians();
        let dist = 100.0;
        let p = [dist * tilt.sin(), 0.0, dist * tilt.cos()];
        let q0 = norm2(&p); // equal conjugates: q = |p|
        let bx = [tilt.cos(), 0.0, -tilt.sin()];
        let by = [0.0, 1.0, 0.0];
        let lambda = 0.0135;

        let n = 10.0;
        let r_guess = 5.0;

        // Solve at theta=0 and theta=PI
        let r_0 = secant_solve(r_guess, 0.0, n, &p, q0, &bx, &by, 0.0, lambda);
        let r_pi = secant_solve(r_guess, std::f64::consts::PI, n, &p, q0, &bx, &by, 0.0, lambda);

        assert!(r_0.is_ok(), "Secant solve failed at theta=0: {:?}", r_0.err());
        assert!(r_pi.is_ok(), "Secant solve failed at theta=PI: {:?}", r_pi.err());

        let r0 = r_0.unwrap();
        let rpi = r_pi.unwrap();

        // Radii should be on the same order of magnitude (within 2x)
        let ratio = r0 / rpi;
        assert!(
            ratio > 0.5 && ratio < 2.0,
            "Radii at theta=0 ({}) and theta=PI ({}) should be comparable, ratio={}",
            r0, rpi, ratio
        );

        // Both should be positive and reasonable
        assert!(r0 > 0.1 && r0 < 50.0, "r(0) = {} out of range", r0);
        assert!(rpi > 0.1 && rpi < 50.0, "r(PI) = {} out of range", rpi);
    }
}
