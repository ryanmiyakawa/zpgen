#!/bin/bash

# Test case for buttressing (gapped zones)
# FSIdx=1 creates gaps between segments (no merging, individual quads)

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
  1 \
  2.0000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-files/buttressing.gds

echo ""
echo "=== Test Complete: Buttressing (Gapped Zones) ==="
echo "This test creates zones with:"
echo "  - FSIdx=1 (gapped zones)"
echo "  - 2um buttress width (gaps between segments)"
echo "  - Each segment exported as individual quad (NO MERGING)"
