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

        // Zone 10 should have radius around 14-15 um for f~100um
        let result = secant_solve(14.0, 0.0, 10.0, &p, q, &bx, &by, 0.0, lambda);

        assert!(result.is_ok());
        let radius = result.unwrap();
        assert!(radius > 10.0 && radius < 20.0, "Radius {} should be between 10 and 20 um", radius);
    }
}
