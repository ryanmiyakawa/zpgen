// Integration test: Generate a real zone plate and verify output
use zpgen_wasm::gds::VecGdsWriter;
use zpgen_wasm::{ZPGenerator, ZPParams};
use std::fs::File;
use std::io::Write;

#[test]
fn test_generate_small_zone_plate() {
    println!("Testing zone plate generation...");

    // Create parameters for a small zone plate
    let params = ZPParams {
        z_tol: 0.01,
        lambda_nm: 13.5,
        p: [0.0, 0.0, 100.0],
        q: 100000.0,
        k_0: [0.0, 0.0, 1.0],
        bx: [1.0, 0.0, 0.0],
        by: [0.0, 1.0, 0.0],
        obscuration_sigma: 0.0,
        na: 0.02, // Small NA for quick test
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

    let mut generator = ZPGenerator::new();

    // Add progress callback
    generator.set_progress_callback(|progress, zone, total| {
        println!("Progress: {:.1}% - Zone {} / {}", progress * 100.0, zone, total);
    });

    let mut writer = VecGdsWriter::new();

    println!("Starting generation...");
    let result = generator.generate(&params, &mut writer);

    match result {
        Ok(_) => {
            println!("✓ Generation successful!");

            let data = writer.into_inner();
            println!("✓ Generated {} bytes of GDS data", data.len());

            // Write to file for inspection
            if let Ok(mut file) = File::create("/tmp/test_zoneplate.gds") {
                file.write_all(&data).unwrap();
                println!("✓ Written to /tmp/test_zoneplate.gds");
            }

            // Verify basic GDS structure
            assert!(data.len() > 110, "GDS file should have header + data");
            assert_eq!(&data[0..2], &[0, 6], "Should start with GDS header");

            println!("✓ All checks passed!");
        }
        Err(e) => {
            panic!("Generation failed: {:?}", e);
        }
    }
}

#[test]
fn test_parameters_validation() {
    println!("Testing parameter variations...");

    // Test with different NA values
    for na in [0.01, 0.02, 0.05] {
        println!("  Testing NA = {}", na);

        let params = ZPParams {
            z_tol: 0.01,
            lambda_nm: 13.5,
            p: [0.0, 0.0, 100.0],
            q: 100000.0,
            k_0: [0.0, 0.0, 1.0],
            bx: [1.0, 0.0, 0.0],
            by: [0.0, 1.0, 0.0],
            obscuration_sigma: 0.0,
            na,
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
        assert!(result.is_ok(), "Should succeed for NA = {}", na);

        let data = writer.into_inner();
        println!("    ✓ Generated {} bytes", data.len());
    }

    println!("✓ All parameter variations passed!");
}
