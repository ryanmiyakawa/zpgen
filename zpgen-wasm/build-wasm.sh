#!/bin/bash
# Build script for WASM target

set -e

echo "Building ZPGen for WebAssembly..."

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "Error: wasm-pack is not installed"
    echo "Install it with: curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh"
    exit 1
fi

# Build for web target
echo "Running wasm-pack build..."
wasm-pack build --target web --release --out-dir examples/browser/pkg

echo "Build complete! Output in examples/browser/pkg/"
echo ""
echo "To test in browser:"
echo "  1. cd examples/browser"
echo "  2. python3 -m http.server 8000"
echo "  3. Open http://localhost:8000 in your browser"
