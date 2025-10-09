#!/bin/bash

# Test case for maximum vertex count (should exceed 200 and trigger buffer flush)
# Uses high NA to create very fine segments, forcing more vertices per zone

../dist/ZPGenHolo \
  0.0100 \
  13.5000 \
  0.000000 0.000000 500.000000 \
  1000000.0000 \
  0.000000 0.000000 1.000000 \
  1.000000 0.000000 0.000000 \
  0.000000 1.000000 0.000000 \
  0.0000 \
  0.2500 \
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
  0.0000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  test-files/max_vertices.gds

echo ""
echo "=== Test Complete: Max Vertices ==="
echo "This test uses high NA (0.25) to create very fine segments"
echo "Check the output for zones that may exceed 200 vertices and require splitting"
