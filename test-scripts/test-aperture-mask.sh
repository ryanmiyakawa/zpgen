#!/bin/bash

# Test case for aperture mask (obscuration mask mode)
# Uses obscuration mask to create complex aperture shapes

../dist/ZPGenHolo \
  0.0100 \
  13.5000 \
  0.000000 0.000000 500.000000 \
  1000000.0000 \
  0.000000 0.000000 1.000000 \
  1.000000 0.000000 0.000000 \
  0.000000 1.000000 0.000000 \
  0.5000 \
  0.0200 \
  0 \
  1 \
  2.0000 \
  45.0000 \
  0.0000 \
  1.0000 \
  0 \
  0.0000 \
  0.0000 \
  10.0000 \
  1 \
  0 \
  0 \
  1 \
  0.0000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-files/aperture_mask.gds

echo ""
echo "=== Test Complete: Aperture Mask ==="
echo "This test creates a zone plate with:"
echo "  - 50% obscuration scale"
echo "  - 45° rotation"
echo "  - 2x anamorphic stretching"
echo "  - Gapped zones (FSIdx=1)"
