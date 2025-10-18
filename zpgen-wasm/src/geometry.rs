// Geometry and pupil masking functions
// Port from zpUtils.cpp

use std::f64::consts::PI;

/// Returns the square of a number
#[inline]
fn square(x: f64) -> f64 {
    x * x
}

/// Checks if the point (cx, cy) lies within the geometric pupil
///
/// Port from zpUtils.cpp lines 146-153
pub fn b_is_in_geometry(cx: f64, cy: f64, obscuration_sigma: f64) -> bool {
    (square(cy) + square(cx) < 1.0)
        && ((square(cy) + square(cx)) >= square(obscuration_sigma))
}

/// Checks if the point (cx, cy) lies within the anamorphic pupil
///
/// Port from zpUtils.cpp lines 139-144
pub fn b_is_in_anamorphic_pupil(
    cx: f64,
    cy: f64,
    anamorphic_fac: f64,
    obscuration_sigma: f64,
) -> bool {
    (square(cy) * square(anamorphic_fac) + square(cx) < 1.0)
        && ((square(cy) * square(anamorphic_fac) + square(cx)) >= square(obscuration_sigma))
}

/// Checks if the point (cx, cy) lies within a custom mask defined by the customMaskIdx
///
/// Port from zpUtils.cpp lines 156-261
#[allow(clippy::excessive_precision)]
pub fn b_is_in_custom_mask(cx: f64, cy: f64, custom_mask_idx: i32) -> bool {
    match custom_mask_idx {
        1 => {
            // Intel MET
            ((cx - 0.65 * (7.5 * PI / 180.0).cos()).powi(2)
                + (cy - 0.65 * (7.5 * PI / 180.0).sin()).powi(2))
            .sqrt()
                <= 0.35
                || ((cx - 0.65 * (127.5 * PI / 180.0).cos()).powi(2)
                    + (cy - 0.65 * (127.5 * PI / 180.0).sin()).powi(2))
                .sqrt()
                    <= 0.35
                || ((cx - 0.65 * (247.5 * PI / 180.0).cos()).powi(2)
                    + (cy - 0.65 * (247.5 * PI / 180.0).sin()).powi(2))
                .sqrt()
                    <= 0.35
        }
        2 => {
            // TDS ZP 2
            let dx = cx * 0.1515;
            let dy = cy * 0.1515;
            square(dx) + square(dy + 0.1515) > square(0.018)
        }
        3 => {
            // TDS ZP 3
            let dx = cx * 0.2396;
            let dy = cy * 0.2396;
            square(dx) / square(0.2396) + square(dy) / square(0.1515) <= 1.0
                && square(dx) + square(dy + 0.1515) > square(0.018)
        }
        4 => {
            // TDS ZP 4
            let dx = cx * 0.285;
            let dy = cy * 0.285;
            square(dx) + square(dy + 0.1515) <= square(0.303)
                && dy + 0.1515 >= 0.0
                && (dy + 0.1515) >= (37.0 * PI / 180.0).tan() * (dx - 0.138)
                && (dy + 0.1515) >= -(37.0 * PI / 180.0).tan() * (dx + 0.138)
                && square(dx) + square(dy + 0.1515) > square(0.018)
        }
        5 => {
            // Square aperture
            let dx = cx / (0.06541 / 0.1);
            let dy = cy / (0.06541 / 0.1);

            (dx.abs() > 0.43 && dx.abs() < 0.86 && dy.abs() > 0.43 && dy.abs() < 0.86)
                || (square(dx) + square(dy) < square(0.2))
        }
        6 => {
            // Diag square aperture
            let dx = cx / (0.06541 / 0.1);
            let dy = cy / (0.06541 / 0.1);

            ((dx + dy).abs() / 2.0_f64.sqrt() > 0.43
                && (dx + dy).abs() / 2.0_f64.sqrt() < 0.86)
                && ((dx - dy).abs() / 2.0_f64.sqrt() > 0.43
                    && (dx - dy).abs() / 2.0_f64.sqrt() < 0.86)
                || (square(dx) + square(dy) < square(0.2))
        }
        7 => {
            // Flip align
            (cx + 0.65).abs() < 0.1 && cy.abs() < 0.25
                || cx.abs() < 0.1 && (cy - 0.65).abs() < 0.25
                || (cx - 0.65).abs() < 0.25 && cy.abs() < 0.1
                || (cx + cy).abs() / 2.0_f64.sqrt() < 0.25
                    && (cx - cy).abs() / 2.0_f64.sqrt() < 0.1
        }
        8 => {
            // Octopole
            ((cx - 0.75 * (PI / 4.0).cos()).powi(2)
                + (cy - 0.75 * (PI / 4.0).sin()).powi(2))
            .sqrt()
                <= 0.1
                || ((cx - 0.75 * (2.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (2.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (3.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (3.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (4.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (4.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (5.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (5.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (6.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (6.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (7.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (7.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || ((cx - 0.75 * (8.0 * PI / 4.0).cos()).powi(2)
                    + (cy - 0.75 * (8.0 * PI / 4.0).sin()).powi(2))
                .sqrt()
                    <= 0.1
                || square(cx) + square(cy) < square(0.2)
        }
        9 => {
            // small rings
            let r = (cx * cx + cy * cy).sqrt();
            let dr = 0.05;
            (r - 0.1).abs() < dr / 2.0
                || (r - 0.3).abs() < dr / 2.0
                || (r - 0.5).abs() < dr / 2.0
                || (r - 0.7).abs() < dr / 2.0
                || (r - 0.9).abs() < dr / 2.0
        }
        10 => {
            // octal rays
            ((cx * 0.0_f64.cos() + cy * 0.0_f64.sin()).abs() < 1.0
                && (cy * 0.0_f64.cos() - cx * 0.0_f64.sin()).abs() < 0.05
                || (cx * (PI / 4.0).cos() + cy * (PI / 4.0).sin()).abs() < 1.0
                    && (-cx * (PI / 4.0).sin() + cy * (PI / 4.0).cos()).abs() < 0.05
                || (cx * (2.0 * PI / 4.0).cos() + cy * (2.0 * PI / 4.0).sin()).abs() < 1.0
                    && (-cx * (2.0 * PI / 4.0).sin() + cy * (2.0 * PI / 4.0).cos()).abs() < 0.05
                || (cx * (3.0 * PI / 4.0).cos() + cy * (3.0 * PI / 4.0).sin()).abs() < 1.0
                    && (-cx * (3.0 * PI / 4.0).sin() + cy * (3.0 * PI / 4.0).cos()).abs()
                        < 0.05)
                && cx * cx + cy * cy > square(0.1)
        }
        11 => {
            // square aperture
            cx.abs() <= 1.0 / 2.0_f64.sqrt() && cy.abs() <= 1.0 / 2.0_f64.sqrt()
        }
        12 => {
            // Horizontal strip:
            cx.abs() <= 1.0
                && cy.abs() <= 0.1
                && !(cx.abs() <= 0.105 && cx.abs() >= 0.085)
        }
        13 => {
            // Vertical strip:
            cy.abs() <= 1.0
                && cx.abs() <= 0.1
                && !(cy.abs() <= 0.105 && cy.abs() >= 0.085)
        }
        15 => {
            // Black ring:
            let r = (cx * cx + cy * cy).sqrt();
            (r <= 1.0) && (r >= 0.105 || r <= 0.085)
        }
        16 => {
            // obscuration only
            false
        }
        17 => {
            // sliver slice
            let r = (cx * cx + cy * cy).sqrt();
            cy.atan2(cx).abs() < 1.0 / PI / 180.0 && r <= 0.99
        }
        18 => {
            // KT iso
            let r = (cx * cx + cy * cy).sqrt();
            r < 1.0 && !((cy + 1.145) * (cy + 1.145) + cx * cx < 0.409 * 0.409)
        }
        _ => true,
    }
}

/// Returns the custom phase for the provided coordinate and mask index
///
/// Port from zpUtils.cpp lines 263-286
pub fn custom_phase(cx: f64, cy: f64, custom_mask_idx: i32) -> f64 {
    match custom_mask_idx {
        1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 => 0.0,
        14 => {
            // spiral phase:
            cy.atan2(cx)
        }
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_b_is_in_geometry() {
        // Point at origin should be in geometry (no obscuration)
        assert!(b_is_in_geometry(0.0, 0.0, 0.0));

        // Point outside unit circle should not be in geometry
        assert!(!b_is_in_geometry(1.5, 0.0, 0.0));

        // Point inside obscuration should not be in geometry
        assert!(!b_is_in_geometry(0.05, 0.0, 0.1));
    }

    #[test]
    fn test_custom_mask_obscuration_only() {
        // Mask 16 should always return false (obscuration only)
        assert!(!b_is_in_custom_mask(0.5, 0.5, 16));
    }

    #[test]
    fn test_custom_phase_spiral() {
        // Mask 14 should return spiral phase
        let phase = custom_phase(1.0, 0.0, 14);
        assert!((phase - 0.0).abs() < 1e-10);

        let phase = custom_phase(0.0, 1.0, 14);
        assert!((phase - PI / 2.0).abs() < 1e-10);
    }
}
