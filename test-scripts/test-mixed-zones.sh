#!/bin/bash

# Test case for mixed zones (FSIdx=2)
# FSIdx=2 creates alternating full zones + gapped zones
# - Odd zones (1,3,5...): Full rings (MERGED into single polygon)
# - Even zones (2,4,6...): Gapped zones (individual quads, NO MERGING)

../dist/ZPGenHolo \
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
  2 \
  2.0000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-files/mixed_zones.gds

echo ""
echo "=== Test Complete: Mixed Zones (FSIdx=2) ==="
echo "This test creates alternating zones:"
echo "  - Odd zones (1,3,5...): Full rings → MERGED polygons"
echo "  - Even zones (2,4,6...): Gapped zones → Individual quads"
echo "  - 2um buttress width in gapped zones"
