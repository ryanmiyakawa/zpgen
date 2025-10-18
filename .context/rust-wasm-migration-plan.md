# Rust + WASM Migration Plan for ZPGenHolo

## Overview

This document outlines the plan to port ZPGenHolo from C++ to Rust compiled as WebAssembly, enabling zone plate generation directly in the browser with streaming output for large GDS files.

## Architecture

```
Rust code → compile to WASM → load in browser → stream GDS chunks → download file
```

### Key Design Principles

1. **Streaming Architecture**: Process zone-by-zone, emit GDS chunks to avoid memory limits
2. **Incremental Migration**: Port in phases with validation at each step
3. **Algorithm Preservation**: Keep same numerical methods and structure as C++ version
4. **Browser Compatibility**: Target 2-4 GB WASM memory limit

## Memory Constraints

### WebAssembly Limits
- Default: 16 MB initial, grows to ~2-4 GB (browser dependent)
- Chrome/Firefox: typically 4 GB max
- Safari: ~2 GB conservative limit

### Expected File Sizes
- Small ZP (NA=0.02, ~100 zones): ~1-5 MB
- Medium ZP (NA=0.05, ~1000 zones): ~20-50 MB
- Large ZP (NA=0.08, ~10k zones): ~100-300 MB

### Solution
The current C++ implementation already processes zone-by-zone (line 373-761 in zpGenHolo.cpp), making it naturally suited for streaming. Port this architecture to Rust with chunk callbacks.

---

## Phase 1: Analysis & Setup

### Task 1.1: Dependency Analysis

**Files to analyze:**
- `src/zpUtils.cpp` and `src/zpUtils.h`
- `src/patternFileUtils.cpp` and `src/patternFileUtils.h`

**Goals:**
- Identify all functions used by zpGenHolo.cpp
- Check for external dependencies beyond C++ stdlib
- Map C++ standard library functions to Rust equivalents
- Document function signatures and purpose

**Expected outcome:**
- Dependency map document
- List of Rust crates needed (likely just `std` and math functions)

### Task 1.2: Streaming Interface Design

**Rust trait for GDS output:**
```rust
/// Trait for streaming GDS file output
pub trait GdsWriter {
    /// Write raw bytes to output
    fn write_chunk(&mut self, data: &[u8]) -> Result<(), GdsError>;

    /// Write GDS header structure
    fn write_header(&mut self, header: &GdsHeader) -> Result<(), GdsError>;

    /// Write a polygon with vertices
    fn write_polygon(&mut self, vertices: &[f64], layer: i32) -> Result<(), GdsError>;

    /// Finalize and close the GDS stream
    fn finalize(&mut self) -> Result<(), GdsError>;
}

/// WASM-specific implementation that calls JavaScript
pub struct StreamingGdsWriter {
    js_callback: js_sys::Function,
    buffer: Vec<u8>,
    buffer_size: usize, // Default: 1MB chunks
}

impl GdsWriter for StreamingGdsWriter {
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
            let array = js_sys::Uint8Array::from(&self.buffer[..]);
            self.js_callback.call1(&JsValue::NULL, &array)
                .map_err(|_| GdsError::JsCallbackFailed)?;
            self.buffer.clear();
        }
        Ok(())
    }
}
```

---

## Phase 2: Rust Port

### Task 2.1: Project Setup

**Create Rust project:**
```bash
cargo new --lib zpgen-wasm
cd zpgen-wasm
```

**Cargo.toml configuration:**
```toml
[package]
name = "zpgen-wasm"
version = "0.1.0"
edition = "2021"

[lib]
crate-type = ["cdylib", "rlib"]

[dependencies]
wasm-bindgen = "0.2"
js-sys = "0.3"
serde = { version = "1.0", features = ["derive"] }
serde-wasm-bindgen = "0.6"

[dev-dependencies]
wasm-bindgen-test = "0.3"

[profile.release]
opt-level = 3
lto = true
```

**Build setup:**
```bash
# Install wasm-pack
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

# Build for web
wasm-pack build --target web --release
```

### Task 2.2: Port Core Numerical Functions

**Files to create:**
- `src/solver.rs` - Port `objectiveFn()` and `secantSolve()`
- `src/optics.rs` - Port optical path length calculations

**Key functions from zpGenHolo.cpp:**

```rust
// Port of objectiveFn (lines 26-46)
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
    // Compute xyz space coordinates from (r, th)
    let u = [r * th.cos(), r * th.sin(), 0.0];

    // Convert to xyz using basis vectors
    let xyz = zp_uxuy_to_xyz(&u, p, bx, by);

    // Compute OPL in waves
    let opd = xyz_to_opl(&xyz, p, q0, lambda) - xyz_to_opl(p, p, q0, lambda);
    let zp_term = -n / 2.0;

    opd + zp_term + phase
}

// Port of secantSolve (lines 48-74)
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
) -> Result<f64, SolverError> {
    const TOL_X: f64 = 0.00001;
    const MAX_ITER: usize = 50;

    let mut r1 = dr_guess;
    let mut r2 = r1 * 1.02;

    for _iter in 0..MAX_ITER {
        let fr1 = objective_fn(r1, th, n, p, q0, bx, by, phase, lambda);
        let fr2 = objective_fn(r2, th, n, p, q0, bx, by, phase, lambda);

        let r0 = r1 - fr1 * (r1 - r2) / (fr1 - fr2);

        // Check convergence
        if (r0 - r1).abs() < TOL_X {
            return Ok(r0);
        }

        // Update guesses
        r2 = r1;
        r1 = r0;
    }

    Err(SolverError::MaxIterationsReached)
}
```

**Validation:**
- Unit tests comparing against known C++ outputs
- Test edge cases (virtual objects, tilted zone plates)

### Task 2.3: Port Coordinate Transformation Utilities

**From zpUtils.cpp:**
- `zpUxUy2XYZ()` - Zone plate coordinates to XYZ
- `freq2zpXYZ()` - Frequency space to XYZ
- `zpRTh2PCCxCy()` - Zone plate polar to pupil coordinates
- `xyz2OPL()` - XYZ to optical path length

**Create `src/transforms.rs`:**
```rust
pub fn zp_uxuy_to_xyz(
    u: &[f64; 3],
    p: &[f64; 3],
    bx: &[f64; 3],
    by: &[f64; 3],
) -> [f64; 3] {
    // Implementation from zpUtils.cpp
}

pub fn xyz_to_opl(
    r: &[f64; 3],
    p: &[f64; 3],
    q0: f64,
    lambda: f64,
) -> f64 {
    // Implementation from zpUtils.cpp
}

// ... other transforms
```

### Task 2.4: Port Geometry and Masking Functions

**From zpUtils.cpp:**
- `bIsInGeometry()` - Check if point is in pupil
- `bIsInCustomMask()` - Custom mask logic (14+ mask types)
- `customPhase()` - Custom phase patterns (spiral, etc.)
- `getPhaseTerm()` - Zernike polynomial evaluation
- `computeZernikePhase()` - Zernike calculations

**Create `src/geometry.rs` and `src/zernike.rs`**

### Task 2.5: Port GDS Export Functions

**From patternFileUtils.cpp:**
- `initGDS()` - Write GDS header
- `exportPolygon()` - Write polygon record
- `renderGDS()` - Write GDS footer

**Create `src/gds.rs`:**
```rust
pub struct GdsHeader {
    pub layer_number: i32,
    pub dbscale: f64,
}

pub struct GdsPolygon {
    pub vertices: Vec<f64>,
    pub layer: i32,
}

impl StreamingGdsWriter {
    pub fn init_gds(&mut self, header: &GdsHeader) -> Result<(), GdsError> {
        // Write GDS header bytes
        // Port from initGDS() in patternFileUtils.cpp
    }

    pub fn export_polygon(&mut self, polygon: &GdsPolygon) -> Result<(), GdsError> {
        // Write polygon record
        // Port from exportPolygon() in patternFileUtils.cpp
    }

    pub fn render_gds(&mut self) -> Result<(), GdsError> {
        // Write GDS footer
        // Port from renderGDS() in patternFileUtils.cpp
        self.flush() // Ensure final chunk is sent
    }
}
```

### Task 2.6: Implement Main Generation Loop

**Create `src/zpgen.rs`:**
```rust
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct ZPGenerator {
    callback: js_sys::Function,
    progress_callback: Option<js_sys::Function>,
}

#[wasm_bindgen]
impl ZPGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(callback: js_sys::Function) -> Self {
        Self {
            callback,
            progress_callback: None,
        }
    }

    #[wasm_bindgen]
    pub fn set_progress_callback(&mut self, callback: js_sys::Function) {
        self.progress_callback = Some(callback);
    }

    #[wasm_bindgen]
    pub fn generate(&mut self, params: JsValue) -> Result<(), JsValue> {
        let params: ZPParams = serde_wasm_bindgen::from_value(params)?;

        let mut writer = StreamingGdsWriter::new(
            self.callback.clone(),
            1024 * 1024, // 1MB buffer
        );

        self.make_zp(params, &mut writer)?;

        Ok(())
    }

    fn make_zp(
        &mut self,
        params: ZPParams,
        writer: &mut impl GdsWriter,
    ) -> Result<(), ZPGenError> {
        // Port main loop from makeZP() (lines 165-848)

        // Initialize GDS
        writer.write_header(&GdsHeader { /* ... */ })?;

        // Compute N_max, N_min (lines 303-344)
        let (n_min, n_max) = self.compute_zone_bounds(&params)?;
        let total_zones = (n_max - n_min) / 2 + 1;

        // Zone generation loop (lines 373-761)
        for n in (n_min..=n_max).step_by(2) {
            // Generate zone n
            self.generate_zone(n, &params, writer)?;

            // Report progress
            if let Some(ref cb) = self.progress_callback {
                let progress = (n - n_min) as f64 / (n_max - n_min) as f64;
                let progress_obj = js_sys::Object::new();
                js_sys::Reflect::set(&progress_obj, &"progress".into(), &progress.into())?;
                js_sys::Reflect::set(&progress_obj, &"zone".into(), &n.into())?;
                js_sys::Reflect::set(&progress_obj, &"totalZones".into(), &total_zones.into())?;
                cb.call1(&JsValue::NULL, &progress_obj)?;
            }
        }

        // Finalize GDS
        writer.finalize()?;

        Ok(())
    }

    fn generate_zone(
        &self,
        n: i32,
        params: &ZPParams,
        writer: &mut impl GdsWriter,
    ) -> Result<(), ZPGenError> {
        // Port zone generation logic (lines 373-761)
        // This is the core algorithm - preserve structure exactly
    }
}

#[derive(serde::Deserialize)]
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
```

---

## Phase 3: JavaScript Integration

### Task 3.1: TypeScript Type Definitions

**Create `types/zpgen.d.ts`:**
```typescript
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
}

export interface ProgressUpdate {
    progress: number; // 0.0 to 1.0
    zone: number;
    totalZones: number;
}

export class ZPGenerator {
    constructor(onChunk: (data: Uint8Array) => void);
    setProgressCallback(onProgress: (update: ProgressUpdate) => void): void;
    generate(params: ZPParams): void;
}

export default function init(module?: WebAssembly.Module): Promise<void>;
```

### Task 3.2: Minimal Working Example

**Create `examples/minimal-example.html`:**
```html
<!DOCTYPE html>
<html>
<head>
    <title>ZPGen WASM Example</title>
</head>
<body>
    <h1>Zone Plate Generator (WASM)</h1>

    <div>
        <label>NA: <input type="number" id="na" value="0.02" step="0.01"></label><br>
        <label>Lambda (nm): <input type="number" id="lambda" value="13.5" step="0.1"></label><br>
        <label>Focal Length (um): <input type="number" id="focal" value="100" step="10"></label><br>
        <button id="generate">Generate Zone Plate</button>
    </div>

    <div id="progress" style="display:none;">
        <progress id="progressBar" max="100" value="0"></progress>
        <span id="progressText"></span>
    </div>

    <script type="module">
        import init, { ZPGenerator } from './pkg/zpgen_wasm.js';

        async function main() {
            await init();

            document.getElementById('generate').addEventListener('click', () => {
                generateZonePlate();
            });
        }

        function generateZonePlate() {
            const na = parseFloat(document.getElementById('na').value);
            const lambda = parseFloat(document.getElementById('lambda').value);
            const focal = parseFloat(document.getElementById('focal').value);

            // Show progress
            document.getElementById('progress').style.display = 'block';

            // Collect GDS chunks
            const chunks = [];

            // Create generator with chunk callback
            const gen = new ZPGenerator((chunk) => {
                chunks.push(chunk);
            });

            // Set progress callback
            gen.setProgressCallback((update) => {
                const percent = Math.round(update.progress * 100);
                document.getElementById('progressBar').value = percent;
                document.getElementById('progressText').textContent =
                    `Zone ${update.zone} / ${update.totalZones} (${percent}%)`;
            });

            // Build parameters
            const params = {
                zTol: 0.01,
                lambda_nm: lambda,
                p: [0, 0, focal],
                q: 100000,
                k_0: [0, 0, 1],
                bx: [1, 0, 0],
                by: [0, 1, 0],
                obscurationSigma: 0,
                na: na,
                nZerns: 0,
                zernikeOrders: [],
                customMaskIdx: 0,
                anamorphicFac: 1,
                anamorphicAzimuth: 0,
                zpcPhase: 0,
                apd: 0,
                apdWindow: 0,
                zpcr2: 0,
                zpcr1: 0,
                biasNm: 10,
                oppositeTone: false,
                randomizeZoneStart: false,
                fsIdx: 0,
                buttressGapWidth: 0,
                buttressPeriod: 0,
                dutyCycle: 0.5,
                layerNumber: 0,
            };

            try {
                // Generate (synchronous in WASM)
                gen.generate(params);

                // Merge chunks and download
                const blob = new Blob(chunks, { type: 'application/octet-stream' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = `zoneplate_na${na}_f${focal}um.gds`;
                a.click();
                URL.revokeObjectURL(url);

                alert('Zone plate generated successfully!');
            } catch (err) {
                alert('Error: ' + err);
            } finally {
                document.getElementById('progress').style.display = 'none';
            }
        }

        main();
    </script>
</body>
</html>
```

### Task 3.3: Web Worker Example (for large files)

**Create `examples/worker-example/zpgen.worker.js`:**
```javascript
import init, { ZPGenerator } from '../../pkg/zpgen_wasm.js';

let initialized = false;

self.onmessage = async (e) => {
    if (e.data.type === 'generate') {
        if (!initialized) {
            await init();
            initialized = true;
        }

        try {
            const chunks = [];

            const gen = new ZPGenerator((chunk) => {
                chunks.push(new Uint8Array(chunk));
                // Send chunk immediately (streaming)
                self.postMessage({
                    type: 'chunk',
                    data: chunk,
                    index: chunks.length - 1,
                });
            });

            gen.setProgressCallback((update) => {
                self.postMessage({
                    type: 'progress',
                    progress: update.progress,
                    zone: update.zone,
                    totalZones: update.totalZones,
                });
            });

            gen.generate(e.data.params);

            self.postMessage({ type: 'complete', totalChunks: chunks.length });
        } catch (err) {
            self.postMessage({ type: 'error', error: err.toString() });
        }
    }
};
```

**Main thread usage:**
```javascript
const worker = new Worker('./zpgen.worker.js', { type: 'module' });

const chunks = [];

worker.onmessage = (e) => {
    switch (e.data.type) {
        case 'chunk':
            chunks[e.data.index] = e.data.data;
            break;

        case 'progress':
            updateProgressBar(e.data.progress, e.data.zone, e.data.totalZones);
            break;

        case 'complete':
            const blob = new Blob(chunks);
            downloadFile(blob, 'zoneplate.gds');
            break;

        case 'error':
            console.error('Worker error:', e.data.error);
            break;
    }
};

// Start generation
worker.postMessage({ type: 'generate', params: zpParams });
```

---

## Phase 4: Testing & Validation

### Task 4.1: Unit Tests

**For each ported module:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secant_solve_simple_case() {
        // Compare against known C++ output
        let result = secant_solve(
            10.0, 0.0, 10.0,
            &[0.0, 0.0, 100.0], 100000.0,
            &[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0],
            0.0, 0.0135
        ).unwrap();

        assert!((result - 14.577).abs() < 0.001);
    }

    #[test]
    fn test_coordinate_transforms() {
        // Test roundtrip: xyz -> pupil -> xyz
    }
}
```

### Task 4.2: Integration Tests

**Test complete generation:**
```rust
#[cfg(test)]
mod integration_tests {
    #[test]
    fn test_small_zp_generation() {
        // Generate small ZP (NA=0.02)
        // Compare binary output with C++ version
        // Should be byte-for-byte identical
    }
}
```

### Task 4.3: Browser Tests

**Use wasm-bindgen-test:**
```rust
#[cfg(target_arch = "wasm32")]
mod browser_tests {
    use wasm_bindgen_test::*;

    #[wasm_bindgen_test]
    fn test_zpgen_in_browser() {
        // Test full generation in browser environment
    }
}
```

### Task 4.4: Performance Validation

**Benchmark targets:**
- Small ZP (NA=0.02, ~100 zones): <1 second
- Medium ZP (NA=0.05, ~1000 zones): <30 seconds
- Large ZP (NA=0.08, ~10k zones): <5 minutes

**Memory targets:**
- Peak memory usage < 500 MB for large ZPs
- Chunk size tuning (test 256KB, 1MB, 4MB)

**Test cases:**
```bash
# Small ZP
./test-scripts/holo-test.sh --na 0.02 --lambda 13.5 --focal 100

# Medium ZP
./test-scripts/holo-test.sh --na 0.05 --lambda 13.5 --focal 5000

# Large ZP
./test-scripts/holo-test.sh --na 0.08 --lambda 13.5 --focal 10000
```

---

## Migration Strategy

### Incremental Approach

1. **Week 1-2: Core Math**
   - Port numerical solvers and transforms
   - No I/O, just pure functions
   - Validate with unit tests

2. **Week 3: GDS Export**
   - Implement streaming GDS writer
   - Test with small generated data

3. **Week 4: Main Loop**
   - Port makeZP() function
   - Generate first complete ZP in Rust
   - Binary comparison with C++

4. **Week 5: WASM Integration**
   - Add wasm-bindgen bindings
   - Create minimal HTML example
   - Test in browser

5. **Week 6: Optimization**
   - Profile performance
   - Optimize hot paths
   - Tune chunk sizes

6. **Week 7: Production**
   - Web Worker integration
   - Error handling
   - User documentation

### Validation Checklist

At each phase:
- [ ] Code compiles without warnings
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Performance meets targets
- [ ] Binary output matches C++ (where applicable)

### Fallback Strategy

Keep C++ version as reference:
- Use for validation
- Fallback if WASM has issues
- Can run C++ via Emscripten as alternative

---

## React Integration Notes

Since you already have a React app, here's how to integrate:

### Install the WASM package

```bash
# In your React app
npm install ./path/to/zpgen-wasm/pkg
```

### Use in React component

```typescript
import { useEffect, useState } from 'react';
import init, { ZPGenerator } from 'zpgen-wasm';

export function ZPGenComponent() {
    const [initialized, setInitialized] = useState(false);
    const [progress, setProgress] = useState(0);
    const [generating, setGenerating] = useState(false);

    useEffect(() => {
        init().then(() => setInitialized(true));
    }, []);

    const handleGenerate = async (params: ZPParams) => {
        setGenerating(true);
        setProgress(0);

        const chunks: Uint8Array[] = [];

        const gen = new ZPGenerator((chunk) => {
            chunks.push(new Uint8Array(chunk));
        });

        gen.setProgressCallback((update) => {
            setProgress(update.progress);
        });

        try {
            gen.generate(params);

            // Download
            const blob = new Blob(chunks);
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'zoneplate.gds';
            a.click();
            URL.revokeObjectURL(url);
        } catch (err) {
            console.error('Generation failed:', err);
        } finally {
            setGenerating(false);
        }
    };

    return (
        <div>
            {/* Your existing UI */}
            {generating && (
                <div>
                    Generating: {Math.round(progress * 100)}%
                </div>
            )}
        </div>
    );
}
```

---

## Next Steps

1. ✅ Read through C++ codebase (completed)
2. ✅ Create migration plan (this document)
3. **TODO**: Analyze zpUtils.cpp and patternFileUtils.cpp dependencies
4. **TODO**: Set up Rust project structure
5. **TODO**: Begin porting core numerical functions

## Resources

- [wasm-bindgen book](https://rustwasm.github.io/wasm-bindgen/)
- [Rust WASM book](https://rustwasm.github.io/book/)
- [Web Workers API](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API)
- Original C++ code: `src/zpGenHolo.cpp`, `src/zpUtils.cpp`, `src/patternFileUtils.cpp`
