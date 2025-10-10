#!/bin/bash
# Test merged polygon refactor with small NA zone plate
# Parameters: zTol lambda_nm px py pz q k0x k0y k0z bxx bxy bxz byx byy byz obscurationSigma NA nZerns customMaskIdx anamorphicFac anamorphicAzimuth ZPCPhase APD APD_window ZPCR2 ZPCR1 bias_nm File_format Opposite_Tone randomizeZoneStart FSIdx buttressGapWidth buttressPeriod block_size NoP IoP blockGrid_pm layerNumber nwaUnitSelection fileName

# Test case: 0.01 NA, 13.5 nm wavelength, no buttresses (FSIdx=0), GDS format
# This should generate merged polygons (one per zone instead of many quads)

../dist/ZPGenHolo \
  0.0100    `# zTol - zone tolerance (lambda/100)` \
  13.5000   `# lambda_nm - wavelength in nm` \
  0.000000 0.000000 500.000000   `# p vector (object position) - 500 um away` \
  1000000.0000  `# q - image distance (1 meter)` \
  0.000000 0.000000 1.000000     `# k_0 vector (illumination direction) - normal incidence` \
  1.000000 0.000000 0.000000     `# bx basis vector - x aligned` \
  0.000000 1.000000 0.000000     `# by basis vector - y aligned` \
  0.0000    `# obscurationSigma - no obscuration` \
  0.0200    `# NA - small numerical aperture for quick test` \
  0         `# nZerns - no Zernike aberrations` \
  0         `# customMaskIdx - no custom mask` \
  1.0000    `# anamorphicFac - no anamorphism` \
  -0.0000   `# anamorphicAzimuth` \
  0.0000    `# ZPCPhase - no phase contrast` \
  1.0000    `# APD - apodization (?)` \
  0         `# APD_window` \
  0.0000    `# ZPCR2` \
  0.0000    `# ZPCR1` \
  10.0000   `# bias_nm - 10 nm bias` \
  1         `# File_format - GDS only` \
  0         `# Opposite_Tone - normal tone` \
  0         `# randomizeZoneStart - start at 0` \
  0         `# FSIdx - NO BUTTRESSES (continuous zones = merged polygons)` \
  0.6000    `# buttressGapWidth` \
  6.0000    `# buttressPeriod` \
  800000    `# block_size` \
  1         `# NoP - number of parts` \
  1         `# IoP - index of part` \
  0         `# blockGrid_pm` \
  1         `# layerNumber` \
  0         `# nwaUnitSelection` \
  test-merged-polygon

echo ""
echo "Test complete! Check output:"
echo "- test-merged-polygon.gds: GDS file with merged polygons"
echo "- test-merged-polygon.txt: Dose information"
echo ""
echo "Expected behavior:"
echo "- Each zone should report '1 merged polygon (X vertices, Y segments)'"
echo "- Vertex count should be ~2*Y (inner + outer edges)"
echo "- Compare file size with old quadrilateral approach"
