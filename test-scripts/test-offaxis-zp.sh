#!/bin/bash

# Test case for off-axis zone plate
# This configuration creates a zone plate for off-axis imaging geometry

./ZPGenHolo \
  0.0100 \
  13.5000 \
  0.000000 0.000000 500.000000 \
  450000.0000 \
  0.104528 0.000000 0.994522 \
  1.000000 0.000000 0.000000 \
  0.000000 1.000000 0.000000 \
  0.0000 \
  0.0825 \
  0 \
  0 \
  1.0000 \
  -0.0000 \
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
  0.6000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-scripts/test-files/offaxis-zp

echo ""
echo "=== Test Complete: Off-Axis Zone Plate ==="
echo "This test creates a zone plate with:"
echo "  - Object distance: 500 μm"
echo "  - Image distance: 450 mm"
echo "  - Off-axis illumination: k0 = [0.104528, 0.000000, 0.994522]"
echo "  - NA: 0.0825"
echo "  - No obscuration"
echo "  - Standard buttressing (W=0.6 dr, T=6 dr)"
