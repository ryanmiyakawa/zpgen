#!/bin/bash

# Test script for large zone plate generation (memory leak testing)
# This generates a large condenser zone plate with high NA to stress test memory usage
# Parameters from Aires condenser configuration - ALL SHAPES IN GEOMETRY

cd "$(dirname "$0")/.."

./src/ZPGenHolo \
    0.1500 \
    13.5000 \
    0.000000 0.000000 100000.000000 \
    100000.0000 \
    0.000000 0.000000 1.000000 \
    1.000000 0.000000 0.000000 \
    0.000000 1.000000 0.000000 \
    0.0000 \
    0.0205 \
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
    1 \
    0 \
    0.6000 \
    6.0000 \
    800000 \
    1 \
    1 \
    0 \
    1 \
    0 \
    Aires_condenser-LR

echo "Large zone plate generation complete"
echo "Check memory usage during execution with Activity Monitor or 'top'"
