// Main zone plate generation logic
// Port from zpGenHolo.cpp

use crate::gds::{export_polygon, init_gds, render_gds, GdsError, GdsWriter};
use crate::geometry::{b_is_in_custom_mask, b_is_in_geometry, custom_phase};
use crate::solver::secant_solve;
use crate::transforms::{cross_product, freq_to_zp_xyz, norm2, xyz_to_opl, zp_rth_to_pc_cxcy};
use crate::zernike::get_phase_term;
use serde::Deserialize;
use std::f64::consts::PI;

#[derive(Debug, Clone, Deserialize)]
pub struct ZPParams {
    pub z_tol: f64,
    pub lambda_nm: f64,
    pub p: [f64; 3],
    pub q: f64,
    pub k_0: [f64; 3],
    pub bx: [f64; 3],
    pub by: [f64; 3],
    pub obscuration_sigma: f64,
    pub na: f64,
    pub n_zerns: usize,
    pub zernike_orders: Vec<f64>,
    pub custom_mask_idx: i32,
    pub anamorphic_fac: f64,
    pub anamorphic_azimuth: f64,
    pub zpc_phase: f64,
    pub apd: f64,
    pub apd_window: f64,
    pub zpcr2: f64,
    pub zpcr1: f64,
    pub bias_nm: f64,
    pub opposite_tone: bool,
    pub randomize_zone_start: bool,
    pub fs_idx: i32,
    pub buttress_gap_width: f64,
    pub buttress_period: f64,
    pub duty_cycle: f64,
    pub layer_number: i32,
}

pub struct ZPGenerator {
    progress_callback: Option<Box<dyn Fn(f64, i32, i32)>>,
}

impl ZPGenerator {
    pub fn new() -> Self {
        Self {
            progress_callback: None,
        }
    }

    pub fn set_progress_callback<F>(&mut self, callback: F)
    where
        F: Fn(f64, i32, i32) + 'static,
    {
        self.progress_callback = Some(Box::new(callback));
    }

    /// Main zone plate generation function
    ///
    /// Port of makeZP() from zpGenHolo.cpp lines 165-848
    pub fn generate<W: GdsWriter>(
        &self,
        params: &ZPParams,
        writer: &mut W,
    ) -> Result<(), GdsError> {
        let lambda = params.lambda_nm / 1000.0; // nm to um
        let bias_um = params.bias_nm / 1000.0; // nm to um

        let mut duty_cycle = params.duty_cycle;
        if duty_cycle == 0.0 {
            duty_cycle = 0.5;
        }

        // Define plane normal
        let n_hat = cross_product(&params.bx, &params.by);

        // Make a mutable copy of p to handle virtual objects
        let mut p = params.p;
        let _virtual_object = if p[2] < 0.0 {
            p[2] = p[2].abs();
            true
        } else {
            false
        };

        // Compute N_max and N_min (maximum and minimum zone numbers)
        let (n_min, n_max) = self.compute_zone_bounds(&params, &p, &n_hat, lambda)?;
        let total_zones = (n_max - n_min) / 2 + 1;

        // Compute focal length for initial radius guesses
        let p_norm = norm2(&p);
        let f = 1.0 / (1.0 / p_norm * p[2].signum() + 1.0 / params.q);

        // Initialize GDS file
        init_gds(writer, params.layer_number)?;

        // Initialize guesses
        let mut r_guess = (n_min as f64 * lambda * f).sqrt();
        let mut r_guess_p1 = ((n_min + 1) as f64 * lambda * f).sqrt();

        let dbscale = 10000.0; // db unit to microns
        const MAX_ZONE_COORDS: usize = 16000;

        // Pre-allocate buffers for merged polygons
        let mut inner_edge_vertices = vec![0.0; MAX_ZONE_COORDS];
        let mut outer_edge_vertices = vec![0.0; MAX_ZONE_COORDS];

        // Zone generation loop
        for n in (n_min..=n_max).step_by(2) {
            let mut use_merged_polygon = false;
            let mut is_first_segment = false;
            let mut inner_vertex_count = 0;
            let mut outer_vertex_count = 0;

            // Compute initial R at an arbitrary angle
            let mut rn = secant_solve(r_guess, 0.0, n as f64, &p, params.q, &params.bx, &params.by, 0.0, lambda)
                .map_err(|e| GdsError::InvalidData(e))?;
            let mut rnp1 = secant_solve(
                r_guess_p1,
                0.0,
                n as f64 + 2.0 * duty_cycle,
                &p,
                params.q,
                &params.bx,
                &params.by,
                0.0,
                lambda,
            )
            .map_err(|e| GdsError::InvalidData(e))?;

            let mut dr = rnp1 - rn;

            // Update guesses for next zone
            r_guess = r_guess_p1;
            r_guess_p1 = r_guess_p1 * (((n + 1) as f64) / (n as f64)).sqrt();

            // Handle buttress modes
            let (is_gap_zone, buttress_width) = match params.fs_idx {
                0 => {
                    // No buttresses
                    if n % 2 == params.opposite_tone as i32 {
                        continue;
                    }
                    (false, 0.0)
                }
                1 => {
                    // Gapped zones
                    if n % 2 == params.opposite_tone as i32 {
                        continue;
                    }
                    (false, dr * (params.buttress_period - params.buttress_gap_width))
                }
                2 => {
                    // Full zones + gaps
                    if n % 2 == params.opposite_tone as i32 {
                        (true, dr * params.buttress_gap_width)
                    } else {
                        (false, 0.0)
                    }
                }
                _ => (false, 0.0),
            };

            // Zone tolerance constraint
            let alpha_zt = 2.0 * (1.0 / (params.z_tol * lambda / rn + 1.0)).acos();

            // Buttress angle
            let alpha_bt = buttress_width / rn;

            // Set alpha to the tighter constraint
            let alpha = if buttress_width != 0.0 && alpha_bt < alpha_zt {
                alpha_bt
            } else if buttress_width == 0.0 {
                alpha_zt
            } else {
                let alpha_ratio = (alpha_bt / alpha_zt) + 1.0;
                alpha_bt * 1.01 / alpha_ratio
            };

            // Determine if this zone uses merged polygons
            use_merged_polygon = !is_gap_zone;

            if use_merged_polygon {
                is_first_segment = true;
            }

            // Start angle (randomization if needed)
            let start_angle = if params.randomize_zone_start {
                // Note: In WASM we'd need a different random source
                0.0
            } else {
                0.0
            };

            let start_angle = if params.custom_mask_idx == 14 {
                start_angle + PI
            } else {
                start_angle
            };

            let mut current_angle = start_angle;
            const ANGLE_EPSILON: f64 = 1e-6;

            // Angle loop
            while current_angle < start_angle + 2.0 * PI + ANGLE_EPSILON {
                // 1) Compute baseline Rn for n (no added phase component)
                rn = secant_solve(rn, current_angle, n as f64, &p, params.q, &params.bx, &params.by, 0.0, lambda)
                    .map_err(|e| GdsError::InvalidData(e))?;
                rnp1 = secant_solve(
                    rnp1,
                    current_angle,
                    n as f64 + 2.0 * duty_cycle,
                    &p,
                    params.q,
                    &params.bx,
                    &params.by,
                    0.0,
                    lambda,
                )
                .map_err(|e| GdsError::InvalidData(e))?;

                // 2) Compute pupil coordinates
                let k_0_2d = [params.k_0[0], params.k_0[1]];
                let pc = zp_rth_to_pc_cxcy(
                    rn,
                    current_angle,
                    &k_0_2d,
                    &p,
                    &params.bx,
                    &params.by,
                    lambda,
                    params.na,
                    params.anamorphic_azimuth,
                );

                let cx = pc[0] * params.anamorphic_fac;
                let cy = pc[1];

                // Compute phase terms
                let mut phase = get_phase_term(
                    cx,
                    cy,
                    &params.zernike_orders,
                    params.n_zerns,
                    params.zpc_phase,
                    params.zpcr1,
                    params.zpcr2,
                );
                phase += custom_phase(cx, cy, params.custom_mask_idx) / (2.0 * PI);

                // Check if point is in geometry
                if !b_is_in_geometry(cx, cy, params.obscuration_sigma)
                    || (params.custom_mask_idx != 0 && !b_is_in_custom_mask(cx, cy, params.custom_mask_idx))
                {
                    // Export accumulated arc if needed
                    if use_merged_polygon && inner_vertex_count > 0 {
                        self.export_accumulated_arc(
                            &inner_edge_vertices[..inner_vertex_count],
                            &outer_edge_vertices[..outer_vertex_count],
                            writer,
                            params.layer_number,
                        )?;
                        inner_vertex_count = 0;
                        outer_vertex_count = 0;
                        is_first_segment = true;
                    }

                    current_angle += alpha;
                    continue;
                }

                // 3) Recompute with phase correction
                if phase != 0.0 {
                    rn = secant_solve(rn, current_angle, n as f64, &p, params.q, &params.bx, &params.by, phase, lambda)
                        .map_err(|e| GdsError::InvalidData(e))?;
                    rnp1 = secant_solve(
                        rnp1,
                        current_angle,
                        n as f64 + 2.0 * duty_cycle,
                        &p,
                        params.q,
                        &params.bx,
                        &params.by,
                        phase,
                        lambda,
                    )
                    .map_err(|e| GdsError::InvalidData(e))?;
                }

                dr = rnp1 - rn;
                let rcm = (rnp1.powi(3) - rn.powi(3)) / (rnp1.powi(2) - rn.powi(2)) * 2.0 / 3.0
                    * (alpha).sin()
                    / alpha;

                // Compute corrected radii at +/- alpha/2
                let rnpa2 = secant_solve(
                    rn,
                    current_angle + alpha / 2.0,
                    n as f64,
                    &p,
                    params.q,
                    &params.bx,
                    &params.by,
                    phase,
                    lambda,
                )
                .map_err(|e| GdsError::InvalidData(e))?;

                let rnma2 = secant_solve(
                    rn,
                    current_angle - alpha / 2.0,
                    n as f64,
                    &p,
                    params.q,
                    &params.bx,
                    &params.by,
                    phase,
                    lambda,
                )
                .map_err(|e| GdsError::InvalidData(e))?;

                // Unbiased CM radii
                let tr1nb = (rcm - dr / 2.0) * 1.0 / (alpha / 2.0).cos();
                let tr2nb = (rcm + dr / 2.0) * 1.0 / (alpha / 2.0).cos();

                // Apply biases
                let (tr1, tr2) = if is_gap_zone {
                    (tr1nb - bias_um / 2.0, tr2nb + bias_um / 2.0)
                } else {
                    (tr1nb + bias_um / 2.0, tr2nb - bias_um / 2.0)
                };

                // Compute corrected radii at the edges
                let tr1ma = tr1 + (rnma2 - rn);
                let tr1pa = tr1 + (rnpa2 - rn);
                let tr2ma = tr2 + (rnma2 - rn);
                let tr2pa = tr2 + (rnpa2 - rn);

                // Trapezoid coordinates in microns
                let trap_coords_um = [
                    tr1ma * (current_angle - alpha / 2.0).cos(),
                    tr1ma * (current_angle - alpha / 2.0).sin(),
                    tr1pa * (current_angle + alpha / 2.0).cos(),
                    tr1pa * (current_angle + alpha / 2.0).sin(),
                    tr2pa * (current_angle + alpha / 2.0).cos(),
                    tr2pa * (current_angle + alpha / 2.0).sin(),
                    tr2ma * (current_angle - alpha / 2.0).cos(),
                    tr2ma * (current_angle - alpha / 2.0).sin(),
                ];

                // Scale to database units
                let mut trap_coords = [0.0; 10];
                for i in 0..8 {
                    trap_coords[i] = dbscale * trap_coords_um[i];
                }
                trap_coords[8] = trap_coords[0]; // Close polygon
                trap_coords[9] = trap_coords[1];

                // Export shape
                if use_merged_polygon {
                    // Accumulate vertices
                    if is_first_segment {
                        inner_edge_vertices[inner_vertex_count] = trap_coords[0];
                        inner_edge_vertices[inner_vertex_count + 1] = trap_coords[1];
                        inner_vertex_count += 2;

                        outer_edge_vertices[outer_vertex_count] = trap_coords[6];
                        outer_edge_vertices[outer_vertex_count + 1] = trap_coords[7];
                        outer_vertex_count += 2;

                        is_first_segment = false;
                    }

                    inner_edge_vertices[inner_vertex_count] = trap_coords[2];
                    inner_edge_vertices[inner_vertex_count + 1] = trap_coords[3];
                    inner_vertex_count += 2;

                    outer_edge_vertices[outer_vertex_count] = trap_coords[4];
                    outer_edge_vertices[outer_vertex_count + 1] = trap_coords[5];
                    outer_vertex_count += 2;

                    // Check vertex limit
                    let total_vertices = (inner_vertex_count / 2) + (outer_vertex_count / 2) + 1;
                    const MAX_GDS_VERTICES: usize = 200;
                    if total_vertices >= MAX_GDS_VERTICES {
                        self.export_accumulated_arc(
                            &inner_edge_vertices[..inner_vertex_count],
                            &outer_edge_vertices[..outer_vertex_count],
                            writer,
                            params.layer_number,
                        )?;
                        inner_vertex_count = 0;
                        outer_vertex_count = 0;
                        is_first_segment = true;
                    }
                } else {
                    // Export individual quad
                    export_polygon(writer, &trap_coords, params.layer_number, 5)?;
                }

                current_angle += alpha;
            }

            // Export final accumulated arc for this zone
            if use_merged_polygon && inner_vertex_count > 0 {
                self.export_accumulated_arc(
                    &inner_edge_vertices[..inner_vertex_count],
                    &outer_edge_vertices[..outer_vertex_count],
                    writer,
                    params.layer_number,
                )?;
            }

            // Report progress
            if let Some(ref callback) = self.progress_callback {
                let progress = (n - n_min) as f64 / (n_max - n_min) as f64;
                callback(progress, n, total_zones);
            }
        }

        // Finalize GDS file
        render_gds(writer)?;

        Ok(())
    }

    /// Compute zone number bounds from NA and sampling
    ///
    /// Port from zpGenHolo.cpp lines 303-344
    fn compute_zone_bounds(
        &self,
        params: &ZPParams,
        p: &[f64; 3],
        n_hat: &[f64; 3],
        lambda: f64,
    ) -> Result<(i32, i32), GdsError> {
        let t_min = lambda / params.na;
        let mut n_max = 1;
        let mut n_min = 1;

        let num_points = 200;
        for i in 0..num_points {
            let th = i as f64 * 2.0 * PI / num_points as f64;
            let fq = [
                1.0 / t_min * th.cos() + params.k_0[0] / lambda,
                1.0 / t_min * th.sin() + params.k_0[1] / lambda,
            ];

            if let Some(r) = freq_to_zp_xyz(&fq, n_hat, p, lambda) {
                let opd = xyz_to_opl(&r, p, params.q, lambda) - xyz_to_opl(p, p, params.q, lambda);

                let temp_n_min = (2.0 * opd) as i32;
                let temp_n_max = (2.0 * opd) as i32 + 1;

                if temp_n_min < n_min {
                    n_min = temp_n_min;
                }
                if temp_n_max > n_max {
                    n_max = temp_n_max;
                }
            }
        }

        Ok((n_min, n_max))
    }

    /// Export accumulated arc polygon
    ///
    /// Port from zpGenHolo.cpp lines 89-160
    fn export_accumulated_arc<W: GdsWriter>(
        &self,
        inner_edge: &[f64],
        outer_edge: &[f64],
        writer: &mut W,
        layer_number: i32,
    ) -> Result<(), GdsError> {
        if inner_edge.is_empty() {
            return Ok(());
        }

        let inner_vertices = inner_edge.len() / 2;
        let outer_vertices = outer_edge.len() / 2;
        let total_vertices = inner_vertices + outer_vertices + 1;

        let mut arc_polygon = vec![0.0; total_vertices * 2];
        let mut arc_idx = 0;

        // Add all inner vertices (CCW)
        for i in 0..inner_edge.len() {
            arc_polygon[arc_idx] = inner_edge[i];
            arc_idx += 1;
        }

        // Add outer vertices in reverse (CW)
        for i in (0..outer_edge.len()).step_by(2).rev() {
            arc_polygon[arc_idx] = outer_edge[i];
            arc_polygon[arc_idx + 1] = outer_edge[i + 1];
            arc_idx += 2;
        }

        // Close polygon
        arc_polygon[arc_idx] = inner_edge[0];
        arc_polygon[arc_idx + 1] = inner_edge[1];

        export_polygon(writer, &arc_polygon, layer_number, total_vertices)?;

        Ok(())
    }
}

impl Default for ZPGenerator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gds::VecGdsWriter;

    #[test]
    fn test_simple_zp_generation() {
        let params = ZPParams {
            z_tol: 0.01,
            lambda_nm: 13.5,
            p: [0.0, 0.0, 100.0],
            q: 100000.0,
            k_0: [0.0, 0.0, 1.0],
            bx: [1.0, 0.0, 0.0],
            by: [0.0, 1.0, 0.0],
            obscuration_sigma: 0.0,
            na: 0.02,
            n_zerns: 0,
            zernike_orders: vec![],
            custom_mask_idx: 0,
            anamorphic_fac: 1.0,
            anamorphic_azimuth: 0.0,
            zpc_phase: 0.0,
            apd: 0.0,
            apd_window: 0.0,
            zpcr2: 0.0,
            zpcr1: 0.0,
            bias_nm: 10.0,
            opposite_tone: false,
            randomize_zone_start: false,
            fs_idx: 0,
            buttress_gap_width: 0.0,
            buttress_period: 0.0,
            duty_cycle: 0.5,
            layer_number: 0,
        };

        let generator = ZPGenerator::new();
        let mut writer = VecGdsWriter::new();

        let result = generator.generate(&params, &mut writer);
        assert!(result.is_ok(), "Generation failed: {:?}", result.err());

        let data = writer.into_inner();
        assert!(data.len() > 102, "GDS file should contain data");
    }
}
