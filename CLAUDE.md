# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ZPGen (Zone Plate Generator) is a C++ application that generates zone plate patterns for nanofabrication. It creates diffractive optical elements by computing zone positions based on optical path length calculations and exports them to various file formats used by electron beam lithography systems.

## Build Commands

### Compile ZPGenHolo (Main Application)
```bash
cd src
clang++ -std=c++11 -g patternFileUtils.cpp zpUtils.cpp zpGenHolo.cpp -o ZPGenHolo
```

Or use the VSCode build task (Cmd+Shift+B) which runs the same command.

### Compile Unit Tests
```bash
cd src
clang++ -std=c++11 -g zpUtils.cpp zpgeomUnitTests.cpp -o zpgeomUnitTests
```

### Run Unit Tests
```bash
cd src
./zpgeomUnitTests
```

## Architecture

### Core Components

**zpGenHolo.cpp** - Main entry point and zone plate generation engine
- `main()` - Parses command-line arguments for all zone plate parameters
- `makeZP()` - Primary zone plate generation function that:
  1. Computes N_max and N_min (maximum/minimum zone numbers) from NA and wavelength
  2. Iterates through each zone number n from N_min to N_max
  3. For each zone, solves for radii using `secantSolve()` at discrete angular positions
  4. Generates trapezoidal segments (or arcs for NWA format)
  5. Exports to selected file format
- `secantSolve()` - Numerical solver that finds zone radii satisfying the optical path condition
- `objectiveFn()` - Objective function representing optical path difference for secant solver

**zpUtils.cpp/h** - Utility functions organized into categories:
- Zernike polynomials: `zgenpt()`, `computeZernikePhase()` for aberration correction
- Geometry functions: `bIsInGeometry()`, `bIsInCustomMask()`, `customPhase()` for pupil masks
- Coordinate transformations between:
  - Spatial frequency (fx, fy) ↔ Zone plate XYZ coordinates
  - Zone plate coordinates ↔ Pupil coordinates (cx, cy)
  - Cartesian XYZ ↔ Zone plate basis vectors (bx, by, bz)
- Optical path computation: `xyz2OPL()` computes path length in wavelengths

**patternFileUtils.cpp/h** - File format export functions
- Supports 5 output formats:
  - Format 0: NWA (ARC) - Heidelberg nanowriter arcs
  - Format 1: GDS - GDSII polygons only
  - Format 2: GDS + TXT - GDSII with dose metadata
  - Format 3: WRV - Pattern file with dose modulation
  - Format 4: GTX - Alternative pattern format
- Each format has `init*()`, `export*()`, and `render*()` functions

### Key Algorithms

**Zone Plate Generation Flow:**
1. Define geometry: p-vector (object position), q (image distance), basis vectors (bx, by), k_0 (illumination direction)
2. Compute zone number bounds (N_min, N_max) by sampling NA circle in frequency space
3. For each zone n:
   - Determine angular step size (alpha) from zone tolerance or buttress constraints
   - Loop through angles theta from 0 to 2π:
     - Solve for radius R_n and R_(n+1) at current angle using secant method
     - Compute pupil coordinates to evaluate masks and aberrations
     - Apply phase corrections from Zernikes and custom masks
     - Re-solve for corrected radii
     - Generate trapezoid coordinates with bias corrections
     - Export shape with dose modulation

**Coordinate System Philosophy:**
- XYZ: Global Cartesian coordinates
- (Ux, Uy): Zone plate coordinates in plane defined by basis vectors bx, by
- (R, theta): Polar coordinates in zone plate plane
- (fx, fy): Spatial frequency coordinates (k-space)
- (cx, cy): Normalized pupil coordinates

### File Organization

- `src/` - Source code and build outputs
  - `zpGenHolo.cpp` - Main zone plate generator
  - `zpUtils.cpp/h` - Mathematical utilities
  - `patternFileUtils.cpp/h` - File format writers
  - `zpgeomUnitTests.cpp` - Unit tests for coordinate transformations
  - `ZPGenHolo` - Compiled executable (gitignored)
  - `*.gds` - Generated output files (gitignored)
- `test-scripts/` - Shell scripts for testing configurations
- `dist/` - Distribution binaries
- `.vscode/` - VSCode build configuration

## Usage Patterns

### Running ZPGenHolo

The executable takes ~30 command-line arguments in a specific order:
```bash
./ZPGenHolo <zTol> <lambda_nm> <px py pz> <q> <k0x k0y k0z> <bxx bxy bxz byx byy byz> \
  <obscurationSigma> <NA> <nZerns> [zernike_orders_and_weights...] <customMaskIdx> \
  <anamorphicFac> <anamorphicAzimuth> <ZPCPhase> <APD> <APD_window> <ZPCR2> <ZPCR1> \
  <bias_nm> <File_format> <Opposite_Tone> <randomizeZoneStart> <FSIdx> \
  <buttressGapWidth> <buttressPeriod> <block_size> <NoP> <IoP> <blockGrid_pm> \
  <layerNumber> <nwaUnitSelection> <fileName>
```

See `test-scripts/holo-test.sh` for example invocations.

### Custom Masks

Custom masks are selected via `customMaskIdx` and defined in `zpUtils.cpp::bIsInCustomMask()`:
- 0: No mask (full aperture)
- 1-13: Various predefined aperture shapes (Intel MET, square, octopole, etc.)
- 14: Spiral phase mask (vortex)
- 15: Ring aperture
- 16: Obscuration only (no zones)
- 17-18: Special test patterns

### Buttress Modes (FSIdx)

- 0: Standard alternating zones (no buttresses)
- 1: Gapped zones - segments within zones
- 2: Full zones with gap zones for mechanical support

## Important Implementation Details

- Zone radii are solved iteratively using secant method to satisfy: OPD(R, theta, n) - OPD(p) = n/2 + phase_correction = 0
- Virtual objects (p[2] < 0) are automatically detected and handled
- Dose modulation is achieved via clockSpeed parameter (0-65535)
- WRV/GTX formats support pixel-aligned gridding via blockGrid_pm parameter
- Zernike polynomials use Fringe/University of Arizona indexing convention
- All spatial units are micrometers internally; wavelength input is nanometers
