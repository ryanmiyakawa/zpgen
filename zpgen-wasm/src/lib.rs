// ZPGen WASM - Zone Plate Generator compiled to WebAssembly
// Port of ZPGenHolo from C++ to Rust

pub mod gds;
pub mod geometry;
pub mod solver;
pub mod transforms;
pub mod zernike;
pub mod zpgen;

// WASM-specific bindings (only compiled for wasm32 target)
#[cfg(target_arch = "wasm32")]
pub mod wasm;

// Re-export main types for convenience
pub use zpgen::{ZPParams, ZPGenerator};

#[cfg(target_arch = "wasm32")]
pub use wasm::WasmZPGenerator;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
