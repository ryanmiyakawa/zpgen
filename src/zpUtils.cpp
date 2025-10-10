 

#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string.h>
#include <vector>
using namespace std;

int nChooseK(int n, int k)
{
    if (k == 0)
    {
        return 1;
    }
    return (n * nChooseK(n - 1, k - 1)) / k;
}

// Retrieves the value the Zernike Polynomial order J at the coordinate [r, th]:
double zgenpt(int j, double r, double th)
{

    // Get dual index [n,m] from j:
    int smct = 2;
    int ct = 0;
    int numInSmct = 3;
    while (j > 0)
    {
        smct += 2;
        j -= numInSmct;
        numInSmct = smct + 1;
        ct++;
    }

    int n = 2 * ct;
    int m = 0;

    for (int q = 1; q <= abs(j); q++)
    {
        if (q % 2 == 1)
        {
            n--;
            m++;
        }
    }
    if ((abs(j)) % 2 == 1)
    {
        m *= -1;
    }

    // Compute value of radial function:
    int p = (n - abs(m)) / 2;
    double kParity = 1;
    double Rv = 0;

    for (int k = 0; k <= p; k += 1)
    {
        if (k % 2 == 1)
            kParity = -1;
        else
            kParity = 1;

        Rv += ((double)(kParity * nChooseK((n - k), k) * nChooseK((n - 2 * k), (p - k)))) * pow(r, (double)(n - 2 * k));
    }

    // Compute value of azimuthal function:
    double Av = 0;
    if (m > 0)
        Av = cos(m * th);
    else if (m < 0)
        Av = -sin(m * th);
    else
        Av = 1;

    return Rv * Av;
}

double square(double x)
{
    return x * x;
}



// sorts coordinates by y value, lowest to highest using bubble sort
void sortRows(double * x, double * y){
    double tempX;
    double tempY;
    for (int m = 0; m < 3; m++){
        for (int k = 0; k < 3; k++){
            if (y[k] > y[k+1]){
                tempY = y[k];
                tempX = x[k];

                y[k] = y[k+1];
                x[k] = x[k+1];

                y[k+1] = tempY;
                x[k+1] = tempX;
            }
        }
    }
}

long roundCoordToPixel(double val, double blockGrid_pixels){

    if (blockGrid_pixels <= 1){
        return (long) round(val);
    }
    // printf("rounding coord: %0.3f -> %0.3f ->")
    return (long) (round( val / blockGrid_pixels ) * blockGrid_pixels);
}


double computeZernikePhase(double cx, double cy, double *orders, int nZerns)
{
    double rn = sqrt(cx * cx + cy * cy);
    double th = atan2(cy, cx);
    double Zv = 0;
    for (int k = 0; k < nZerns; k++)
    {
        Zv += orders[nZerns + k] * zgenpt((int)orders[k], rn, th);
    }
    return Zv;
}

double getPhaseTerm(double cx, double cy, double *orders, int nZerns, double ZPCPhase, double ZPCR1, double ZPCR2)
{
    double r = sqrt(cx * cx + cy * cy);
    double ph = 0;
    if (r <= ZPCR1 && r >= ZPCR2)
        ph = ph + ZPCPhase / (2 * M_PI);

    ph = ph + computeZernikePhase(cx, cy, orders, nZerns);
    return ph;
}

bool bIsInAnamorphicPupil(double cx, double cy, double anamorphicFac, double obscurationSigma)
{

    return (square(cy) * square(anamorphicFac) + square(cx) < 1) &&
           ((square(cy) * square(anamorphicFac) + square(cx)) >= square(obscurationSigma));
}

bool bIsInGeometry(double cx, double cy, double obscurationSigma)
{
    // if ((square(cy)  + square(cx)) >= square(obscurationSigma)){
    //     printf("cx: %0.3f, cy: %0.3f, rad: %0.5f\n", cx, cy, sqrt(square(cy)  + square(cx)));
    // }
     return (square(cy)  + square(cx) < 1) &&
           ((square(cy)  + square(cx)) >= square(obscurationSigma));
}

// Specify custom mask inputs
bool bIsInCustomMask(double cx, double cy, int customMaskIdx)
{
    double dx, dy, r, dr;
    switch (customMaskIdx){
        case 1: // Intel MET
            return (sqrt(square(cx - .65*cos(7.5*M_PI/180)) + square(cy - .65*sin(7.5*M_PI/180))) <= .35) ||
                        (sqrt(square(cx - .65*cos(127.5*M_PI/180)) + square(cy - .65*sin(127.5*M_PI/180))) <= .35) ||
                        (sqrt(square(cx - .65*cos(247.5*M_PI/180)) + square(cy - .65*sin(247.5*M_PI/180))) <= .35);
            
        case 2: // TDS ZP 2
            dx = cx * 0.1515;
            dy = cy * 0.1515;
            return square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            
        case 3: // TDS ZP 3
            dx = cx * 0.2396;
            dy = cy * 0.2396;
            return square(dx)/square(0.2396) + square(dy)/square(0.1515) <= 1 &&  // outer ellipse
                    square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            
        case 4: // TDS ZP 4
            dx = cx * 0.285;
            dy = cy * 0.285;
            return square(dx) + square(dy + 0.1515) <= square(0.303) &&  // outer parent
                    dy + 0.1515 >= 0 &&  // x-axis
                    (dy + 0.1515) >= tan(37*M_PI/180)*(dx - 0.138) && // diagonal lines
                    (dy + 0.1515) >= -tan(37*M_PI/180)*(dx + 0.138) &&
                    square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            
        case 5: // Square aperture
            dx = cx / (.06541/.1);
            dy = cy / (.06541/.1);
            
            return (abs(dx) > 0.43 && abs(dx) < 0.86 &&
                    abs(dy) > 0.43 && abs(dy) < 0.86)||
            (square(dx) + square(dy) < square(0.2));
            
        case 6: // Diag square aperture
            dx = cx / (.06541/.1);
            dy = cy / (.06541/.1);
            
            return ((abs(dx + dy)/sqrt(2) > 0.43 && abs(dx + dy)/sqrt(2) < 0.86) &&
            (abs(dx - dy)/sqrt(2) > 0.43 && abs(dx - dy)/sqrt(2) < 0.86)) ||
            (square(dx) + square(dy) < square(0.2));
            
        case 7: // Flip align
            return   (abs(cx + .65) < .1 && abs(cy) < .25) ||
            (abs(cx) < .1 && abs(cy - .65) < .25) ||
            (abs(cx - 0.65) < .25 && abs(cy) < .1) ||
            (abs(cx + cy)/sqrt(2) < .25 && abs(cx - cy)/sqrt(2) < .1);
       
        case 8: // Octopole
            return  (sqrt(square(cx - 0.75*cos(M_PI/4))   + square(cy - 0.75*sin(M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(2*M_PI/4)) + square(cy - 0.75*sin(2*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(3*M_PI/4)) + square(cy - 0.75*sin(3*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(4*M_PI/4)) + square(cy - 0.75*sin(4*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(5*M_PI/4)) + square(cy - 0.75*sin(5*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(6*M_PI/4)) + square(cy - 0.75*sin(6*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(7*M_PI/4)) + square(cy - 0.75*sin(7*M_PI/4))) <= 0.1) ||
                    (sqrt(square(cx - 0.75*cos(8*M_PI/4)) + square(cy - 0.75*sin(8*M_PI/4))) <= 0.1) ||
                    square(cx) + square(cy) < square(0.2);
            
        case 9: // small rings
            r = sqrt(cx*cx + cy*cy);
            dr = 0.05;
            return (abs(r - .1) < dr/2 ||
                    abs(r - .3) < dr/2 ||
                    abs(r - .5) < dr/2 ||
                    abs(r - .7) < dr/2 ||
                    abs(r - .9) < dr/2);
            break;
        case 10: // octal rays
            return ((abs(cx*cos(0) + cy*sin(0)) < 1 && abs(cy*cos(0) - cx*sin(0)) < .05) ||
                    (abs(cx*cos(1*M_PI_4) + cy*sin(1*M_PI_4)) < 1 && abs(-cx*sin(1*M_PI_4) + cy*cos(1*M_PI_4)) < .05) ||
                    (abs(cx*cos(2*M_PI_4) + cy*sin(2*M_PI_4)) < 1 && abs(-cx*sin(2*M_PI_4) + cy*cos(2*M_PI_4)) < .05) ||
                    (abs(cx*cos(3*M_PI_4) + cy*sin(3*M_PI_4)) < 1 && abs(-cx*sin(3*M_PI_4) + cy*cos(3*M_PI_4)) < .05)) &&
                    cx*cx + cy*cy > square(.1);
        case 11: // square aperture
            return abs(cx) <= 1/sqrt(2) && abs(cy) <= 1/sqrt(2);
            break;
        case 12: // Horizontal strip:
            return abs(cx) <= 1 && abs(cy) <= 0.1 && !(abs(cx) <= 0.105 && abs(cx) >= 0.085);
            break;
        case 13: // Vertical strip:
            return abs(cy) <= 1 && abs(cx) <= 0.1 && !(abs(cy) <= 0.105 && abs(cy) >= 0.085);
            break;
        case 15: // Black ring:
            r = sqrt(cx*cx + cy*cy);
            return (r <= 1) && (r >= 0.105 || r <= 0.085);
            break;

        case 16: // obscuration only
            return 0;
        case 17: // sliver slice
            r = sqrt(cx*cx + cy*cy);            
            return  abs(atan2(cy, cx)) < M_1_PI/180;// && (r <= 0.99);
            return 0;
        case 18: // KT iso 
            r = sqrt(cx*cx + cy*cy);    

            return r < 1 && !((cy + 1.145)*(cy + 1.145) + cx*cx < 0.409*0.409);   
            break;

    }
    return true;
}

double customPhase(double cx, double cy, int customMaskIdx)
{
    switch (customMaskIdx)
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
        return 0;

    case 14: // spiral phase:
        return atan2(cy, cx);
    }
    return 0;
}


void getZPCoordinatesFromFrequency(double kx, double ky, double P, double *coordinates, double beta)
{
    // Define k-vector:
    double kz = sqrt(1 - kx * kx - ky * ky);

    double p0x = 0;
    double p0y = 0;
    double p0z = P;

    double n0x = 0;
    double n0y = sin(beta);
    double n0z = cos(beta);

    double d = (p0x * n0x + p0y * n0y + p0z * n0z) / (kx * n0x + ky * n0y + kz * n0z);

    coordinates[0] = d * kx;
    coordinates[1] = d * ky / cos(beta);
}

void getPupilCoordinatesFromZPCoordinates(double x, double yp, double P, double *coordinates, double beta, double NA, double CRA, double anamorphicFac)
{
    // Transform [xp, yp, zp] to x-z rectilinear coords:

    double y = yp * cos(beta);
    double z = P - yp * sin(beta);

    // [x,y,z] is then our k-vector, so normalize it:
    double absK = sqrt(x * x + y * y + z * z);

    double kx = x / absK;
    double ky = y / absK;

    coordinates[0] = kx / NA;
    coordinates[1] = anamorphicFac * (ky - sin(CRA)) / NA;
}

/**
 * a2, b2, and c2 are the squared sides
 */
double computeTriangleAreaByHeron(double a2, double b2, double c2){

    return sqrt( 4*a2 * b2 - (a2 + b2 - c2) * (a2 + b2 - c2)) / 4;
}

double computeAverageRadiusDiff(double * trapCoords, double dbScale, double blockGrid_pixels, 
                                double * drs, long * coords){
    double x1, y1, x2, y2, x3, y3, x4, y4, r1, r2;

    // first triangle:
    x1 = ((double) roundCoordToPixel(trapCoords[0], blockGrid_pixels))/dbScale;
    y1 = ((double) roundCoordToPixel(trapCoords[1], blockGrid_pixels))/dbScale;
    x2 = ((double) roundCoordToPixel(trapCoords[2], blockGrid_pixels))/dbScale;
    y2 = ((double) roundCoordToPixel(trapCoords[3], blockGrid_pixels))/dbScale;
    x3 = ((double) roundCoordToPixel(trapCoords[4], blockGrid_pixels))/dbScale;
    y3 = ((double) roundCoordToPixel(trapCoords[5], blockGrid_pixels))/dbScale;
    x4 = ((double) roundCoordToPixel(trapCoords[6], blockGrid_pixels))/dbScale;
    y4 = ((double) roundCoordToPixel(trapCoords[7], blockGrid_pixels))/dbScale;

    r1 = (sqrt(square(x1) + square(y1)) + sqrt(square(x2) + square(y2)))/2 ;
    r2 = (sqrt(square(x3) + square(y3)) + sqrt(square(x4) + square(y4)))/2 ;

    drs[0] = r2;
    drs[1] = r1;

    for (int k = 0; k < 8; k++){
        coords[k] = roundCoordToPixel(trapCoords[k], blockGrid_pixels);
    }

    return r2 - r1;
}


double computeAreaOfOrderedQuadrilateral(double * trapCoords, double dbScale, double blockGrid_pixels){

    double x1, y1, x2, y2, x3, y3, x4, y4, A1, A2;

    // first triangle:
    x1 = ((double) roundCoordToPixel(trapCoords[0] * dbScale, blockGrid_pixels))/dbScale;
    y1 = ((double) roundCoordToPixel(trapCoords[1] * dbScale, blockGrid_pixels))/dbScale;
    x2 = ((double) roundCoordToPixel(trapCoords[2] * dbScale, blockGrid_pixels))/dbScale;
    y2 = ((double) roundCoordToPixel(trapCoords[3] * dbScale, blockGrid_pixels))/dbScale;
    x3 = ((double) roundCoordToPixel(trapCoords[4] * dbScale, blockGrid_pixels))/dbScale;
    y3 = ((double) roundCoordToPixel(trapCoords[5] * dbScale, blockGrid_pixels))/dbScale;

//  printf("1: [%0.2f], 1: [ %0.2f], 1: [ %0.2f]\n", square(x2 - x1) + square(y2 - y1), square(x3 - x2) + square(y3 - y2), square(x1 - x3) + square(y1 - y3));


    A1 = computeTriangleAreaByHeron( 
            square(x2 - x1) + square(y2 - y1),
            square(x3 - x2) + square(y3 - y2),
            square(x1 - x3) + square(y1 - y3)
        );
    // printf("A1: %0.10f\n", A1);
    // second triangle:
    x3 = ((double) roundCoordToPixel(trapCoords[4] * dbScale, blockGrid_pixels))/dbScale;
    y3 = ((double) roundCoordToPixel(trapCoords[5] * dbScale, blockGrid_pixels))/dbScale;
    x4 = ((double) roundCoordToPixel(trapCoords[6] * dbScale, blockGrid_pixels))/dbScale;
    y4 = ((double) roundCoordToPixel(trapCoords[7] * dbScale, blockGrid_pixels))/dbScale;
    x1 = ((double) roundCoordToPixel(trapCoords[0] * dbScale, blockGrid_pixels))/dbScale;
    y1 = ((double) roundCoordToPixel(trapCoords[1] * dbScale, blockGrid_pixels))/dbScale;

    A2 = computeTriangleAreaByHeron( 
            square(x4 - x3) + square(y4 - y3),
            square(x1 - x4) + square(y1 - y4),
            square(x3 - x1) + square(y3 - y1)
        );

    // printf("[%0.4f,%0.4f,%0.4f,%0.4f; %0.4f,%0.4f,%0.4f,%0.4f] \n", 
    //     x1, x2, x3, x4, y1, y2, y3, y4);


    return (A1 + A2) ;
}

long getPolyClockSpeedFromRadiusRatio(double minAreaFraction, double maxAreaFraction, double dr, double * trapCoords,
    double dbScale, double blockGrid_pixels ){

    double drs[2];
    long coords[8];
    double griddedDr = computeAverageRadiusDiff(trapCoords, dbScale, blockGrid_pixels, drs, coords);

    double areaFraction = dr / griddedDr;


    if (areaFraction < minAreaFraction){
            areaFraction = areaFraction;
    }
    if (areaFraction > maxAreaFraction){
        areaFraction = maxAreaFraction;
    }

    long clockSpeed = (long) (65535 - ((areaFraction - minAreaFraction)/(maxAreaFraction - minAreaFraction) * 65535));

//    printf("dr rat: %0.6f/%0.6f, drs: [%0.6f, %0.6f], coords: [%ld, %ld;%ld, %ld;%ld, %ld;%ld, %ld] clock speed: %ld\n",
//      dr, griddedDr, drs[0], drs[1], coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], coords[6], coords[7], clockSpeed);


    // printf("dr rat: %0.6f/%0.6f, drs: [%0.6f, %0.6f], coords: [%0.3f, %0.3f;%0.3f, %0.3f;%0.3f, %0.3f;%0.3f, %0.3f]clock speed: %ld\n",
    //  dr, griddedDr, drs[0], drs[1], coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], coords[6], coords[7], clockSpeed);

    return clockSpeed;

}

long getPolyClockSpeedFromAreaFraction(double minClockSpeed, double maxClockSpeed, double * trapCoords, double * trapCoordsNoBias,
    double dbScale, double blockGrid_pixels ){

    double griddedArea = computeAreaOfOrderedQuadrilateral(trapCoords, dbScale, blockGrid_pixels);
    double nominalArea = computeAreaOfOrderedQuadrilateral(trapCoordsNoBias, dbScale, 0);


    double thisClockSpeed = griddedArea / nominalArea;

    // printf("Area fraction %0.3f = %0.3f/%0.3f.  [m,M] = [%0.3f, %0.3f]\n", thisClockSpeed, griddedArea*1000, nominalArea*1000, minClockSpeed, maxClockSpeed);

    if (thisClockSpeed < minClockSpeed){
            thisClockSpeed = minClockSpeed;
    }
    if (thisClockSpeed > maxClockSpeed){
        thisClockSpeed = maxClockSpeed;
    }

    long clockSpeed = (long) ((thisClockSpeed - minClockSpeed)/(maxClockSpeed - minClockSpeed) * 65535);

    return clockSpeed;
}

// Computes requisite clock speed to perform dose bias
long getPolyClockSpeed(double dR1, double dRN, double dRn, double bias)
{
    double maxClockSpeed = (dRN - bias) / dRN;
    double minClockSpeed = (dR1 - bias) / dR1;

    // True when bias = 0, so set all doses to maximum dose
    if (maxClockSpeed == minClockSpeed)
    {
        return (long)65535;
    }

    double thisClockSpeed = (dRn - bias) / dRn;

    // Require thisClockSpeed to be within bounds (inequalities are backward due to definition)
    if (thisClockSpeed > minClockSpeed)
    {
        thisClockSpeed = minClockSpeed;
    }
    if (thisClockSpeed < maxClockSpeed)
    {
        thisClockSpeed = maxClockSpeed;
    }

    long clockSpeed = (long)((maxClockSpeed - thisClockSpeed) / (maxClockSpeed - minClockSpeed) * 65535);
    

    return clockSpeed;
}



// ZPGEOM Utils:
double norm2(double * v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double scalarDotProduct(double *v1, double *v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void normVector(double * v){
    double norm = norm2(v);
    v[0] = v[0]/norm;
    v[1] = v[1]/norm;
    v[2] = v[2]/norm;
}

// Compute the cross product of two vectors
void crossProduct( double a[3],  double b[3], double result[3]) {
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}


static inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double norm3(const double v[3]) {
    return sqrt(dot3(v, v));
}


// Coverts spatial frequency defiend by (fx,fy) to a point in xyz space
bool freq2zpXYZ(const double f[2],
                           const double n_in[3],
                           const double p_on_plane[3],
                           double lambda,
                           double out[3]) {

    const double fx = f[0];
    const double fy = f[1];

    // 1) Build unit direction l from spatial frequencies
    const double s2 = (lambda*fx)*(lambda*fx) + (lambda*fy)*(lambda*fy);
    if (s2 > 1.0) {
        // Evanescent: no real propagation direction
        return false;
    }
    double lz = sqrt(fmax(0.0, 1.0 - s2));
    double l[3] = { lambda*fx, lambda*fy, lz }; // already unit length

    // 2) Normalize plane normal (without mutating input)
    double n[3] = { n_in[0], n_in[1], n_in[2] };
    const double nmag = norm3(n);
    if (nmag == 0.0) return false;
    n[0] /= nmag; n[1] /= nmag; n[2] /= nmag;

    // 3) Choose lz sign so intersection is in front of the origin (t >= 0) if possible
    // Compute denom with current sign
    double denom = dot3(l, n);
    double numer = dot3(p_on_plane, n);

    if (denom == 0.0) {
        // Try flipping lz once; if still zero, it's parallel
        l[2] = -l[2];
        denom = dot3(l, n);
        if (denom == 0.0) return false;
    }

    double t = numer / denom;

    // If you require forward intersection only (t >= 0), try flipping lz once to see if that helps
    if (t < 0.0) {
        l[2] = -l[2];
        denom = dot3(l, n);
        if (denom == 0.0) return false;
        t = numer / denom;
        if (t < 0.0) {
            // Plane is behind the chosen propagation direction
            return false;
        }
    }

    // 4) Intersection point
    out[0] = t * l[0];
    out[1] = t * l[1];
    out[2] = t * l[2];
    return true;


    // double fx = f[0];
    // double fy = f[1];

    // // Make sure u and n are unit vectors:
    // double p0 = norm2(p);
    // double u[3] = {p[0]/p0, p[1]/p0, p[2]/p0};
    // normVector(n);

    // // build l-hat vector:
    // double l_norm = 1/lambda * sqrt(1-((fx * lambda)*(fx * lambda) + (fy * lambda)*(fy * lambda)));
    // double l[3] = {fx, fy, l_norm};
    // normVector(l);

    // // flip lz if necessary due to ambiguity of square root:
    // if (p[2] < 0){
    //     l[2] = -l[2];
    // }

    // // solve p0 = l_0 + d*l for d:
    // double d = p0 * scalarDotProduct(u, n) / scalarDotProduct(l, n);

    // result[0] = d*l[0];
    // result[1] = d*l[1];
    // result[2] = d*l[2];
}

void zpCoord2KVector(double r[3], double k[2]) {
    double r_norm = norm2(r);
    k[0] = r[0] / r_norm ;
    k[1] = r[1] / r_norm ;
}

void zpXYZ2UxUy(const double r[3], const double p[3],  double bx[3],  double by[3], double result[3]) {
    // Normalization of basis vectors is skipped because assuming they are normalized
    
    // Point in Cartesian coordinates wrt p0
    double rp[3] = {r[0] - p[0], r[1] - p[1], r[2] - p[2]};

    // compute bz from bx and by
    double bz[3];
    crossProduct(bx, by, bz);

    // Project onto basis
    result[0] = scalarDotProduct(rp, bx);
    result[1] = scalarDotProduct(rp, by);
    result[2] = scalarDotProduct(rp, bz);
}

void zpUxUy2XYZ( double U[3],  double p0[3],  double bx[3],  double by[3], double * result) {
    // Normalization of basis vectors is skipped because assuming they are normalized

    // Convert U to Cartesian and write in terms of p0 origin
    result[0] = p0[0] + U[0]*bx[0] + U[1]*by[0];
    result[1] = p0[1] + U[0]*bx[1] + U[1]*by[1];
    result[2] = p0[2] + U[0]*bx[2] + U[1]*by[2];
}

/**
 * @brief computes the pathlength in waves from the origin to r, then to the point q.  Handles
 * Virtual and Real sources and objects
 * 
 * @param r_o - point in xyz space on lens coming from object
 * @param p - vector that points along optical axis (from object to image)
 * @param q0 - distance along optical axis from zp to image
 * @return double 
 */
double xyz2OPL( double r_o[3], double p[3], double q0, double lambda){

    double p_norm = norm2(p);
    double p_hat[3] = { p[0]/p_norm, p[1]/p_norm, p[2]/p_norm };

    // distance from object (origin) to point r_o
    double L1 = norm2(r_o);

    // projection of r_o onto the optical axis
    double proj = r_o[0]*p_hat[0] + r_o[1]*p_hat[1] + r_o[2]*p_hat[2];

    // how much farther along p we need to go to reach the image plane
    double delta = q0 - proj;

    // compute image point r_i = r_o + delta * p_hat
    double r_i[3] = {
        r_o[0] + delta * p_hat[0],
        r_o[1] + delta * p_hat[1],
        r_o[2] + delta * p_hat[2]
    };

    // distance from r_o to r_i
    double L2_vec[3] = {
        r_i[0] - r_o[0],
        r_i[1] - r_o[1],
        r_i[2] - r_o[2]
    };
    double L2 = norm2(L2_vec);

    // return total optical path length in waves
    return (L1 + L2) / lambda;
    
    
    // double p_norm = norm2(p);

    // // vector from image to point xyz on zp.
    // double * r_i = new double[3];

    // // vector relationship between r_o and r_i: r_i = p + q - r_o
    // // q and p point in the same direction: q = q0 * p/|p|
    // for (int i = 0; i < 3; i++)
    //     r_i[i] = p[i] + q0 * p[i] / p_norm - r_o[i];
    

    // return (norm2(r_o) + norm2(r_i)) / lambda;
}

// Composite function useful for directly grabbing pupil coordinates
// Takes a point defined in ZP coordinates and computes pupil coordinates
void zpRTh2PCCxCy(double R, double th, double * k_0, double * p, double * bx, double * by, double lambda, double NA, double * C, double anamorphicAzimuth){
    double ux = R * cos(th);
    double uy = R * sin(th);
    double * k = new double[2]; 
    double * U = new double[2];
    double * Cxy = new double[2];
    U[0] = ux;
    U[1] = uy;

    double * r = new double[3];
    zpUxUy2XYZ(U, p, bx, by, r);
    zpCoord2KVector(r, k);

    Cxy[0] = (k[0] - k_0[0]) / NA;
    Cxy[1] = (k[1] - k_0[1]) / NA;

    // Determine in-plane angle of k by taking atan2(k[1],k[0]):
    double k_azi = anamorphicAzimuth; // -atan2(k_0[1], k_0[0]);

    // Rotate C by k_azi:
    C[0] = cos(k_azi) * Cxy[0] - sin(k_azi) * Cxy[1];
    C[1] = sin(k_azi) * Cxy[0] + cos(k_azi) * Cxy[1];

}

// Converts spatial frequency in space to a point defined by [ux, uy] in ZP coordinates using basis vectors {bx, by}
void freq2zpUxUy(double * fq, double * n, double * p, double * bx, double * by, double lambda, double * U){
    
    // F -> r
    double * r = new double[3];
    freq2zpXYZ(fq, n, p, lambda, r);

    // r -> U
    zpXYZ2UxUy(r, p, bx, by, U);
}
