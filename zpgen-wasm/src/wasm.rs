// WASM bindings for browser usage
// Provides JavaScript-accessible interface to ZPGen

use crate::gds::{GdsError, GdsWriter};
use crate::zpgen::{ZPGenerator, ZPParams};
use wasm_bindgen::prelude::*;

/// WASM-specific GDS writer that calls back to JavaScript
pub struct WasmGdsWriter {
    callback: js_sys::Function,
    buffer: Vec<u8>,
    buffer_size: usize,
}

impl WasmGdsWriter {
    pub fn new(callback: js_sys::Function, buffer_size: usize) -> Self {
        Self {
            callback,
            buffer: Vec::with_capacity(buffer_size),
            buffer_size,
        }
    }
}

impl GdsWriter for WasmGdsWriter {
    fn write_chunk(&mut self, data: &[u8]) -> Result<(), GdsError> {
        self.buffer.extend_from_slice(data);

        // Flush buffer when it reaches threshold
        if self.buffer.len() >= self.buffer_size {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), GdsError> {
        if !self.buffer.is_empty() {
            // Create Uint8Array from buffer
            let array = js_sys::Uint8Array::new_with_length(self.buffer.len() as u32);
            array.copy_from(&self.buffer);

            // Call JavaScript callback
            self.callback
                .call1(&JsValue::NULL, &array)
                .map_err(|_| GdsError::InvalidData("JavaScript callback failed".to_string()))?;

            self.buffer.clear();
        }
        Ok(())
    }
}

/// WASM-exposed zone plate generator
#[wasm_bindgen]
pub struct WasmZPGenerator {
    inner: ZPGenerator,
    chunk_callback: Option<js_sys::Function>,
    progress_callback: Option<js_sys::Function>,
}

#[wasm_bindgen]
impl WasmZPGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(chunk_callback: js_sys::Function) -> Self {
        Self {
            inner: ZPGenerator::new(),
            chunk_callback: Some(chunk_callback),
            progress_callback: None,
        }
    }

    #[wasm_bindgen(js_name = setProgressCallback)]
    pub fn set_progress_callback(&mut self, callback: js_sys::Function) {
        self.progress_callback = Some(callback);
    }

    /// Compute zone information (count, bounds, parent NA) without generating geometry
    /// This is a lightweight calculation for UI preview
    #[wasm_bindgen(js_name = computeZoneCount)]
    pub fn compute_zone_count(&self, params_js: JsValue) -> Result<JsValue, JsValue> {
        // Deserialize parameters from JavaScript
        let params: ZPParams = serde_wasm_bindgen::from_value(params_js)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse parameters: {:?}", e)))?;

        // Compute zone info
        let zone_info = self.inner
            .compute_zone_count(&params)
            .map_err(|e| JsValue::from_str(&format!("Zone count calculation failed: {:?}", e)))?;

        // Serialize to JavaScript object
        serde_wasm_bindgen::to_value(&zone_info)
            .map_err(|e| JsValue::from_str(&format!("Failed to serialize zone info: {:?}", e)))
    }

    #[wasm_bindgen]
    pub fn generate(&mut self, params_js: JsValue) -> Result<(), JsValue> {
        // Deserialize parameters from JavaScript
        let params: ZPParams = serde_wasm_bindgen::from_value(params_js)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse parameters: {:?}", e)))?;

        // Log deserialized params so we can verify what WASM actually received
        web_sys::console::log_1(&format!(
            "[WASM] generate() params: p={:?}, q={}, na={}, lambda_nm={}, duty_cycle={}",
            params.p, params.q, params.na, params.lambda_nm, params.duty_cycle
        ).into());

        // Set up progress callback
        if let Some(ref cb) = self.progress_callback {
            let cb_clone = cb.clone();
            self.inner.set_progress_callback(move |progress, zone, total| {
                let obj = js_sys::Object::new();
                js_sys::Reflect::set(&obj, &"progress".into(), &progress.into()).unwrap();
                js_sys::Reflect::set(&obj, &"zone".into(), &zone.into()).unwrap();
                js_sys::Reflect::set(&obj, &"totalZones".into(), &total.into()).unwrap();
                let _ = cb_clone.call1(&JsValue::NULL, &obj);
            });
        }

        // Create WASM writer
        let chunk_cb = self.chunk_callback.clone().ok_or_else(|| {
            JsValue::from_str("No chunk callback set")
        })?;

        let mut writer = WasmGdsWriter::new(chunk_cb, 1024 * 1024); // 1MB chunks

        // Generate zone plate
        self.inner
            .generate(&params, &mut writer)
            .map_err(|e| JsValue::from_str(&format!("{:?}", e)))?;

        Ok(())
    }
}

// Export types for TypeScript
#[wasm_bindgen(typescript_custom_section)]
const TS_APPEND_CONTENT: &'static str = r#"
export interface ZPParams {
    zTol: number;
    lambda_nm: number;
    p: [number, number, number];
    q: number;
    k_0: [number, number, number];
    bx: [number, number, number];
    by: [number, number, number];
    obscurationSigma: number;
    na: number;
    nZerns: number;
    zernikeOrders: number[];
    customMaskIdx: number;
    anamorphicFac: number;
    anamorphicAzimuth: number;
    zpcPhase: number;
    apd: number;
    apdWindow: number;
    zpcr2: number;
    zpcr1: number;
    biasNm: number;
    oppositeTone: boolean;
    randomizeZoneStart: boolean;
    fsIdx: number;
    buttressGapWidth: number;
    buttressPeriod: number;
    dutyCycle: number;
    layerNumber: number;
    maxGdsVertices: number;
}

export interface ProgressUpdate {
    progress: number;
    zone: number;
    totalZones: number;
}
"#;
