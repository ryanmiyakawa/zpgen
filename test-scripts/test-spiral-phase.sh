#!/bin/bash

# Test case for spiral phase mask (customMaskIdx=14)
# This tests the fix for merged polygons with geometry discontinuities

./ZPGenHolo \
  0.0100 \
  13.5000 \
  0.000000 0.000000 500.000000 \
  1000000.0000 \
  0.000000 0.000000 1.000000 \
  1.000000 0.000000 0.000000 \
  0.000000 1.000000 0.000000 \
  0.0000 \
  0.0200 \
  0 \
  14 \
  1.0000 \
  0.0000 \
  0.0000 \
  1.0000 \
  0 \
  0.0000 \
  0.0000 \
  10.0000 \
  1 \
  0 \
  0 \
  0 \
  0.0000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-scripts/test-files/spiral_phase

echo ""
echo "=== Test Complete: Spiral Phase ==="
echo "This test creates a zone plate with:"
echo "  - Spiral phase mask (customMaskIdx=14)"
echo "  - No obscuration"
echo "  - Tests arc export with phase discontinuity"
