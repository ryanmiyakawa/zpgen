#ifndef ZP_UTILS_H
#define ZP_UTILS_H

#include <cmath>
#include <vector>

// n Choose k function
int nChooseK(int n, int k);

// Retrieves the value of the Zernike Polynomial order J at the coordinate [r, th]
double zgenpt(int j, double r, double th);

// Returns the square of a number
double square(double x);

// sorts coordinates by y value, lowest to highest using bubble sort
void sortRows(double * x, double * y);

// Rounds the coordinate value to the nearest pixel
long roundCoordToPixel(double val, double blockGrid_pixels);

// Computes the phase of Zernike polynomial
double computeZernikePhase(double cx, double cy, double *orders, int nZerns);

// Returns phase term for the provided parameters
double getPhaseTerm(double cx, double cy, double *orders, int nZerns, double ZPCPhase, double ZPCR1, double ZPCR2);

// Checks if the point (cx, cy) lies within the anamorphic pupil
bool bIsInAnamorphicPupil(double cx, double cy, double anamorphicFac, double obscurationSigma);

// Checks if the point (cx, cy) lies within the geometric shape
bool bIsInGeometry(double cx, double cy, double obscurationSigma);

// Checks if the point (cx, cy) lies within a custom mask defined by the customMaskIdx
bool bIsInCustomMask(double cx, double cy, int customMaskIdx);

// Returns the custom phase for the provided coordinate and mask index
double customPhase(double cx, double cy, int customMaskIdx);

void getZPCoordinatesFromFrequency(double kx, double ky, double P, double *coordinates, double beta);

// Computes and returns the ZPC coordinates from the provided frequency
void getPupilCoordinatesFromZPCoordinates(double x, double yp, double P, double *coordinates, double beta, double NA, double CRA, double anamorphicFac);

double computeTriangleAreaByHeron(double a2, double b2, double c2);

double computeAverageRadiusDiff(double * trapCoords, double dbScale, double blockGrid_pixels, 
                                double * drs, long * coords);

double computeAreaOfOrderedQuadrilateral(double * trapCoords, double dbScale, double blockGrid_pixels);

long getPolyClockSpeedFromRadiusRatio(double minAreaFraction, double maxAreaFraction, double dr, double * trapCoords,    double dbScale, double blockGrid_pixels );

long getPolyClockSpeedFromAreaFraction(double minClockSpeed, double maxClockSpeed, double * trapCoords, double * trapCoordsNoBias,
    double dbScale, double blockGrid_pixels );

long getPolyClockSpeed(double dR1, double dRN, double dRn, double bias);


// Utility functions for ZPGEOM
double norm2(double * v);
double scalarDotProduct(double *v1, double *v2);
void normVector(double * v);
void crossProduct( double a[3],  double b[3], double result[3]);

// Main functions
bool freq2zpXYZ(const double f[2],
                           const double n_in[3],
                           const double p_on_plane[3],
                           double lambda,
                           double out[3]);
void zpCoord2KVector(double r[3], double result[2]);
void zpXYZ2UxUy(const double r[3], const double p[3],  double bx[3],  double by[3],  double bz[3], double result[3]);
void zpUxUy2XYZ( double U[3],  double p0[3],  double bx[3],  double by[3], double result[3]);
double xyz2OPL( double r[3], double p[3], double q0, double lambda);

// Composite function useful for directly grabbing pupil coordinates
void zpRTh2PCCxCy(double R, double th, double * k, double * p, double * bx, double * by, double lambda, double NA, double * C, double anamorphicAzimuth);
void freq2zpUxUy(double * fq, double * n, double * p, double * bx, double * by, double lambda, double * U);

#endif // ZP_UTILS_H
