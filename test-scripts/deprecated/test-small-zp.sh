#!/bin/bash
# Test script based on verified working parameters
# 0.02 NA, 13.5 nm wavelength, 500 um object distance, 1000 mm image distance

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
  0 \
  0.6000 \
  6.0000 \
  800000 \
  1 \
  1 \
  0 \
  1 \
  0 \
  small-zp

echo ""
echo "Test complete! Output files:"
echo "- small-zp.gds"
echo "- small-zp.txt (if format was GDS+txt)"
