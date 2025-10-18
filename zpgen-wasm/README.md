# ZPGen WASM - Zone Plate Generator for WebAssembly

Rust port of ZPGenHolo (C++ zone plate generator) compiled to WebAssembly for browser-based zone plate generation.

## Project Status

✅ **Phase 1: Core Rust Port (COMPLETE)**
- All C++ code has been ported to Rust
- Modules implemented:
  - `solver.rs` - Numerical secant solver for zone radii
  - `transforms.rs` - Coordinate transformations (XYZ ↔ pupil ↔ frequency)
  - `zernike.rs` - Zernike polynomial evaluation
  - `geometry.rs` - Pupil masks and custom geometries (18 mask types)
  - `gds.rs` - GDS file format export with streaming support
  - `zpgen.rs` - Main zone plate generation loop

⏳ **Next Steps:**
- Add WASM bindings with wasm-bindgen
- Create browser example
- Add streaming callback for large files
- Performance testing and optimization

## Architecture

### Streaming Design

The Rust implementation maintains the C++ zone-by-zone processing architecture, making it ideal for streaming output:

```rust
pub trait GdsWriter {
    fn write_chunk(&mut self, data: &[u8]) -> Result<(), GdsError>;
    fn flush(&mut self) -> Result<(), GdsError>;
}
```

For WASM, we'll implement:
```rust
pub struct StreamingGdsWriter {
    js_callback: js_sys::Function,
    buffer: Vec<u8>,
    buffer_size: usize,  // 1MB chunks
}
```

### Key Differences from C++

1. **Memory Safety**: Rust's ownership system eliminates pointer arithmetic and manual memory management
2. **Error Handling**: Uses `Result<T, E>` instead of C-style error codes
3. **Streaming**: Trait-based design allows multiple output backends (Vec, WASM callback, file)
4. **Testing**: Built-in unit tests in each module

## Module Overview

### solver.rs
Port of `objectiveFn()` and `secantSolve()` from zpGenHolo.cpp:
- Numerical root-finding for zone radii
- Tolerance: 0.00001 μm, max 50 iterations

### transforms.rs
Port of zpUtils.cpp coordinate transformation functions:
- `zpuxuy_to_xyz()` - Zone plate coordinates → XYZ
- `xyz_to_opl()` - Optical path length calculation
- `freq_to_zp_xyz()` - Spatial frequency → XYZ
- `zp_rth_to_pc_cxcy()` - Zone plate polar → pupil coordinates

### zernike.rs
Port of Zernike polynomial evaluation:
- Fringe/University of Arizona indexing convention
- Supports arbitrary orders for aberration correction

### geometry.rs
Port of pupil masking and custom geometry functions:
- 18 custom mask types (Intel MET, octopole, spiral phase, etc.)
- Anamorphic pupil support
- Obscuration handling

### gds.rs
GDS binary file format writer:
- Header/footer generation
- Polygon export with layer support
- Streaming-ready trait-based design

### zpgen.rs
Main zone plate generation logic:
- Zone number bounds calculation from NA
- Zone-by-zone polygon generation
- Merged polygon support (reduces vertex count by ~1000x)
- Progress callbacks
- Buttress modes (free-standing, gapped zones)

## Build Instructions

```bash
# Standard Rust build
cargo build --release

# Run tests
cargo test

# Build for WASM (requires wasm-pack)
wasm-pack build --target web --release
```

## Testing

Each module includes unit tests:

```bash
# Run all tests
cargo test

# Run specific module tests
cargo test --test solver
cargo test --test transforms
```

## Usage Example (Native Rust)

```rust
use zpgen_wasm::{ZPGenerator, ZPParams};
use zpgen_wasm::gds::VecGdsWriter;

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
    // ... other parameters
};

let generator = ZPGenerator::new();
let mut writer = VecGdsWriter::new();

generator.generate(&params, &mut writer)?;

let gds_data = writer.into_inner();
// Write to file or send to browser
```

## Dependencies

- `serde` - Parameter serialization for WASM
- `wasm-bindgen` - Rust ↔ JavaScript interop
- `js-sys` - JavaScript standard library bindings

## Performance Targets

Based on C++ benchmarks:
- Small ZP (NA=0.02, ~100 zones): <1 second
- Medium ZP (NA=0.05, ~1000 zones): <30 seconds
- Large ZP (NA=0.08, ~10k zones): <5 minutes

Memory usage: <500 MB peak for large zone plates

## License

Copyright © 2017-2025 Ryan Miyakawa and Henry Wang. All rights reserved.
