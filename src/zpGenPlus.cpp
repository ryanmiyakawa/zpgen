/**
 * 
 *  zpGenPlus.cpp
 *  ZPGenPlus
 *  Created by Ryan Miyakawa on 2/17/17.
 *  Copyright © 2017 Ryan Miyakawa and Henry Wang. All rights reserved.
 * 
 * 
 * 
 */


#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string.h>

using namespace std;

int nChooseK(int n, int k){
    if (k == 0) {
        return 1;
    }
    return (n * nChooseK(n - 1, k - 1)) / k;
}

// Retrieves the value the Zernike Polynomial order J at the coordinate [r, th]:
double zgenpt(int j, double r, double th) {
    
    // Get dual index [n,m] from j:
    int smct = 2;
    int ct = 0;
    int numInSmct = 3;
    while (j > 0) {
        smct += 2;
        j -= numInSmct;
        numInSmct = smct + 1;
        ct++;
    }
    
    int n = 2 * ct;
    int m = 0;
    
    for (int q = 1; q <= abs(j); q++) {
        if (q % 2 == 1) {
            n--;
            m++;
        }
    }
    if ((abs(j)) % 2 == 1) {
        m *= -1;
    }
    
    // Compute value of radial function:
    int p = (n - abs(m)) / 2;
    double kParity = 1;
    double Rv = 0;
    
    for (int k = 0; k <= p; k += 1) {
        if (k % 2 == 1)
            kParity = -1;
        else
            kParity = 1;
        
        Rv += ((double)(kParity*nChooseK((n - k), k)
            *nChooseK((n - 2 * k), (p - k))))
            *pow(r, (double)(n - 2 * k));
    }

    // Compute value of azimuthal function:
    double Av = 0;
    if (m > 0)
        Av = cos(m*th);
    else if (m < 0)
        Av = -sin(m*th);
    else
        Av = 1;
    
    return Rv*Av;
}

double square(double x){
    return x*x;
}

void encode32(long aCoord, char * cPart){
    cPart[0] = (aCoord >> 24) & 255;
    cPart[1] = (aCoord >> 16) & 255;
    cPart[2] = (aCoord >> 8) & 255;
    cPart[3] = (aCoord) & 255;
}

void encodePoly32(long * coords, char * cCoords){
    char cPart[4];
    for (int k = 0; k < 10; k++){
        encode32(coords[k], cPart);
        for (int m = 0; m < 4; m++){
            cCoords[k*4 + m] = cPart[m];
        }
    }
}

// sorts coordinates by y value, lowest to highest using bubble sort
void sortRows(float * x, float * y){
    float tempX;
    float tempY;
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
        return (long) roundf(val);
    }
    return (long) (roundf( val/blockGrid_pixels ) * blockGrid_pixels);
}


// Fractures trapeizoids into two triangles and a horizontal trapezoid for WRV
void fractureAndWriteWRVPoly(long * coords, FILE * outputFile, long clockSpeed, double blockGrid_pixels){
    float heightTol = 5; // shapes with smaller than 0.5 nm height are ignored
    // Need to sort coords
    float x[4] = {(float)coords[0], (float)coords[2], (float)coords[4], (float)coords[6]};
    float y[4] = {(float)coords[1], (float)coords[3], (float)coords[5], (float)coords[7]};
    
    // printf("Presorting: ");
    // printf("[[%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]]\n", x[0], y[0],x[1], y[1],x[2], y[2],x[3], y[3]);

    // Sortrows
    sortRows(x,y);
    // printf("Postsorting: ");
    // printf("[[%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]]", x[0], y[0],x[1], y[1],x[2], y[2],x[3], y[3]);
    
    // Find intersections.  Two cases: y[3] -> y[1] (more common), or y[3] -> y[0] (less common)
    float xU, xD;

    // There are two types of trapezoids:
    /**
     *  1) 0 -> 3 is max slope ore min slope
     *  2) 0 -> 3 is in between 0->2 and 0->1
     */

    float m03 = atan2((y[3] - y[0]),(x[3] - x[0]));
    float m02 = atan2((y[2] - y[0]),(x[2] - x[0]));
    float m01 = atan2((y[1] - y[0]),(x[1] - x[0]));

    if ( (m03 > m01) == (m02 > m03)){  
        // printf(": Used type 1 trapezoid (mid) \n");
        xU = x[3] - (x[3] - x[1])*(y[3] - y[2])/(y[3] - y[1]);
        xD = x[2] - (x[2] - x[0])*(y[2] - y[1])/(y[2] - y[0]);
    } else { // case 2 trapezoid
    //  printf(": Used type 1 trapezoid (L/R) \n");
        // xU = x[3] - (x[3] - x[0])*(y[3] - y[2])/(y[3] - y[0]);
        // xD = x[3] - (x[3] - x[0])*(y[3] - y[1])/(y[3] - y[0]);
        xU = x[0] + (x[3] - x[0])*(y[2] - y[0])/(y[3] - y[0]);
        xD = x[0] + (x[3] - x[0])*(y[1] - y[0])/(y[3] - y[0]);
    }
    float xDL, xDR, xUL, xUR;
    
    // Find orientation of lower and upper triangles
    if (xD < x[1]){
        xDL = xD;
        xDR = x[1];
    } else {
        xDL = x[1];
        xDR = xD;
    }
    if (xU < x[2]){
        xUL = xU;
        xUR = x[2];
    } else {
        xUL = x[2];
        xUR = xU;
    }
    
//     printf("Lower tri: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
//             (long) x[0], (long)y[0], (long)xDR, (long)y[1], (long)x[0], (long)xDL); // lower tri
//  printf("Upper tri: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
//             (long) xUL, (long)y[2], (long)x[3], (long)y[3], (long)xUR, (long)x[3]); // upper tri
// printf("Middle trap: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
//                 (long) xDL, (long)y[1], (long)xUR, (long)y[2], (long)xDR, (long)xUL); // middle trap

    if (y[1] - y[0] > heightTol){
        fprintf(outputFile, "Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
            roundCoordToPixel(x[0], blockGrid_pixels), 
            roundCoordToPixel(y[0], blockGrid_pixels), 
            roundCoordToPixel(xDR, blockGrid_pixels), 
            roundCoordToPixel(y[1], blockGrid_pixels), 
            roundCoordToPixel(x[0], blockGrid_pixels), 
            roundCoordToPixel(xDL, blockGrid_pixels)
        ); // lower tri
    }

    if (y[3] - y[2] > heightTol){
        fprintf(outputFile, "Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
            roundCoordToPixel(xUL, blockGrid_pixels), 
            roundCoordToPixel(y[2], blockGrid_pixels),
            roundCoordToPixel(x[3], blockGrid_pixels), 
            roundCoordToPixel(y[3], blockGrid_pixels), 
            roundCoordToPixel(xUR, blockGrid_pixels),  
            roundCoordToPixel(x[3], blockGrid_pixels)
        ); // upper tri
    }
    if (y[2] - y[1] > heightTol){
        fprintf(outputFile, "Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
            roundCoordToPixel(xDL, blockGrid_pixels), 
            roundCoordToPixel(y[1], blockGrid_pixels), 
            roundCoordToPixel(xUR, blockGrid_pixels), 
            roundCoordToPixel(y[2], blockGrid_pixels), 
            roundCoordToPixel(xDR, blockGrid_pixels), 
            roundCoordToPixel(xUL, blockGrid_pixels)
        ); // middle trap
    }
}

void exportArc(double R, double dR, double theta, double dTheta, double dose, double nwaUnit,
               double offsetX, double offsetY, FILE * outputFile){
    fprintf(outputFile, "ARC_S %12.4f %10.4f %10.4f %10.4f  %5.2f  %5.2f %5.2f %5.2f D%3.4f S1.00 T-14000\n",
            R/nwaUnit, fmod((theta * 180/M_PI), 360), dR/nwaUnit, dTheta * 180/M_PI,
            offsetX/nwaUnit, offsetY/nwaUnit, 0., 0., dose);

}


void exportPolygon(long * coords, unsigned char * polyPre, unsigned char * polyPost,
                   unsigned char * polyForm, FILE * outputFile, int File_format, long clockSpeed, double blockGrid_pixels){
    switch (File_format){
        case 1:
        case 2:
            char cCoords[40];
            fwrite(polyPre, sizeof(char), 16, outputFile);
            fwrite(polyForm, sizeof(char), 4, outputFile);
            encodePoly32(coords, cCoords);
            fwrite(cCoords, sizeof(char), 40, outputFile);
            fwrite(polyPost, sizeof(char), 4, outputFile);
            break;
        case 3:
            fractureAndWriteWRVPoly(coords, outputFile, clockSpeed, blockGrid_pixels);
            break;
        }
}

double computeZernikePhase(double cx, double cy, double * orders, int nZerns){
    double rn = sqrt(cx*cx + cy*cy);
    double th = atan2(cy, cx);
    double Zv = 0;
    for (int k = 0; k < nZerns; k++) {
        Zv += orders[nZerns + k] * zgenpt((int)orders[k], rn, th);
    }
    return Zv;
}

double getPhaseTerm(double cx, double cy, double * orders, int nZerns, double ZPCPhase, double ZPCR1, double ZPCR2){
    double r = sqrt(cx*cx + cy*cy);
    double ph = 0;
    if (r <= ZPCR1 && r >= ZPCR2)
        ph = ph + ZPCPhase/(2*M_PI);
    
    ph = ph + computeZernikePhase(cx, cy, orders, nZerns);
    return ph;
}

bool bIsInGeometry(double cx, double cy, double obscurationSigma, double CRA){
    double r = sqrt(cx*cx + cy*cy);
    if (CRA == 0){
        return (r >= obscurationSigma);
    }else {
        return (r >= obscurationSigma && r <= 1.002);
    }

}

// Specify custom mask inputs
bool bIsInCustomMask(double cx, double cy, int customMaskIdx){
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

    }
    return true;
}

double customPhase(double cx, double cy, int customMaskIdx){
    switch (customMaskIdx){
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

bool bIsInAnamorphicPupil(double cx, double cy, double anamorphicFac, double obscurationSigma){
    
    return (square(cy)*square(anamorphicFac) + square(cx) < 1) &&
    ((square(cy)*square(anamorphicFac) + square(cx)) >= square(obscurationSigma));
}

double objectiveFn(double r, double th, double N, double p, double q,
                   double phase, double lambda, double beta, int virtualObject) {
    
    double pp, qp, rp, zpTerm, phTerm, plTermP, plTermQ;
    
    pp = p - r*sin(beta)*sin(th);
    qp = q + r*sin(beta)*sin(th);
    rp = r*sqrt(sin(th)*sin(th)*cos(beta)*cos(beta) + cos(th) * cos(th));
    
    plTermP = rp*sqrt(pp*pp/rp/rp + 1) - pp;
    plTermQ = rp*sqrt(qp*qp/rp/rp + 1) - qp;
    
    // Zero out possible precision errors
    if (plTermP < 1e-8){
        plTermP = 0;
    }
    if (plTermQ < 1e-8){
        plTermQ = 0;
    }
    zpTerm = -N*lambda/2;
    phTerm = phase*lambda;
    
    return plTermP + virtualObject*plTermQ + zpTerm + phTerm;
}

double secantSolve(double dRGuess, double th, double N, double p, double q, double phase, double lambda, double beta, int virtualObject){
    double tolX = 0.00001;
    int maxIter = 50;
    
    //Make guesses:
    double x1, x2, fxm1, fxm2, R0;
    x1 = dRGuess;
    x2 = x1*1.02;
    for (int currentIter = 0; currentIter < maxIter; currentIter++){
        fxm1 = objectiveFn(x1, th, N, p, q, phase, lambda, beta, virtualObject);
        fxm2 = objectiveFn(x2, th, N, p, q, phase, lambda, beta, virtualObject);

        R0 = x1 - fxm1*(x1 - x2) / (fxm1 - fxm2);
        if (abs(R0 - x1) < tolX){
            // printf("p: %0.6f, q: %0.6f, N: %d, lambda: %0.3f, beta: %0.2f, ph: %0.3f, isVirt: %d\n", p, q, (int)N, lambda, beta, phase, virtualObject);
            return R0;

        }
    
        //Set new guesses
        x2 = x1;
        x1 = R0;
    }
    printf("MAXIMUM ITERATIONS REACHED\n");
    return x1;
}


// Computes dose ratio from zone width bias
double getDoseFromBias(double dr, double bias){
    return dr/(dr - bias);
}

/**
 * 
 * 
 * 
 */
void getZPCoordinatesFromFrequency(double kx, double ky, double P, double * coordinates, double beta){
    // Define k-vector:
    double kz = sqrt(1 - kx*kx - ky*ky);

    double p0x = 0;
    double p0y = 0;
    double p0z = P;

    double n0x = 0;
    double n0y = sin(beta);
    double n0z = cos(beta);

    double d = (p0x*n0x + p0y*n0y + p0z*n0z)/(kx*n0x + ky*n0y + kz*n0z);

    coordinates[0] = d * kx;
    coordinates[1] = d * ky / cos(beta);
}

void getPupilCoordinatesFromZPCoordinates(double x, double yp, double P, double * coordinates, double beta, double NA, double CRA){
    // Transform [xp, yp, zp] to x-z rectilinear coords:

    double y = yp*cos(beta);
    double z = P - yp*sin(beta);

    // [x,y,z] is then our k-vector, so normalize it:
    double absK = sqrt(x*x + y*y + z*z);

    double kx = x/absK;
    double ky = y/absK;

    coordinates[0] = kx/NA;
    coordinates[1] = (ky - sin(CRA))/NA;
}

/**
 * a2, b2, and c2 are the squared sides
 */
double computeTriangleAreaByHeron(double a2, double b2, double c2){

    return sqrt( 4*a2 * b2 - (a2 + b2 - c2) * (a2 + b2 - c2)) / 4;
}



double computeAreaOfOrderedQuadrilateral(double * trapCoords, double dbScale, double blockGrid_pixels){

    double x1, y1, x2, y2, x3, y3, A1, A2;

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
    x1 = ((double) roundCoordToPixel(trapCoords[4] * dbScale, blockGrid_pixels))/dbScale;
    y1 = ((double) roundCoordToPixel(trapCoords[5] * dbScale, blockGrid_pixels))/dbScale;
    x2 = ((double) roundCoordToPixel(trapCoords[6] * dbScale, blockGrid_pixels))/dbScale;
    y2 = ((double) roundCoordToPixel(trapCoords[7] * dbScale, blockGrid_pixels))/dbScale;
    x3 = ((double) roundCoordToPixel(trapCoords[0] * dbScale, blockGrid_pixels))/dbScale;
    y3 = ((double) roundCoordToPixel(trapCoords[1] * dbScale, blockGrid_pixels))/dbScale;

    A2 = computeTriangleAreaByHeron( 
            square(x2 - x1) + square(y2 - y1),
            square(x3 - x2) + square(y3 - y2),
            square(x1 - x3) + square(y1 - y3)
        );

    return (A1 + A2) ;
}

long getPolyClockSpeedFromAreaFraction(double maxClockSpeed, double minClockSpeed, double * trapCoords, double nominalArea,
    double dbScale, double blockGrid_pixels ){


    double thisClockSpeed = computeAreaOfOrderedQuadrilateral(trapCoords, dbScale, blockGrid_pixels) / nominalArea;

    // printf("Area fraction %0.2f = %0.3f/%0.3f\n", thisClockSpeed, computeAreaOfOrderedQuadrilateral(trapCoords) , nominalArea);
    if (thisClockSpeed > minClockSpeed){
            thisClockSpeed = minClockSpeed;
    }
    if (thisClockSpeed < maxClockSpeed){
        thisClockSpeed = maxClockSpeed;
    }

    long clockSpeed = (long) ((maxClockSpeed - thisClockSpeed)/(maxClockSpeed - minClockSpeed) * 65535);

    return clockSpeed;
}

// Computes requisite clock speed to perform dose bias
long getPolyClockSpeed(double dR1, double dRN, double dRn, double bias){
    double maxClockSpeed = (dRN - bias)/dRN;
    double minClockSpeed = (dR1 - bias)/dR1;
    
    // True when bias = 0, so set all doses to maximum dose
    if (maxClockSpeed == minClockSpeed){
        return (long) 65535;
    }
    
    double thisClockSpeed = (dRn - bias)/dRn;

    // Require thisClockSpeed to be within bounds (inequalities are backward due to definition)
    if (thisClockSpeed > minClockSpeed){
        thisClockSpeed = minClockSpeed;
    }
    if (thisClockSpeed < maxClockSpeed){
        thisClockSpeed = maxClockSpeed;
    }


    long clockSpeed = (long) ((maxClockSpeed - thisClockSpeed)/(maxClockSpeed - minClockSpeed) * 65535);
    

    return clockSpeed;
}

void writeSupportTxtHeader(double dR1, double dRN, double bias, FILE * supportTextFile){
    fprintf(supportTextFile, "Dose %f %f\n", getDoseFromBias(dRN, bias), getDoseFromBias(dR1, bias));
}

void writeSupportTxtZoneDose(long clockSpeed, int n, double Rn,  FILE * supportTextFile){
    fprintf(supportTextFile, "Zone %d %ld %f\n", n, clockSpeed, Rn);
}

void initWRV(FILE * outputFile, double minDose, double maxDose, long block_size, long block_size_unit_pm){
    fprintf(outputFile, "patdef %ld %ld %ld %f %f 0 0\n", block_size_unit_pm, block_size, block_size, maxDose, minDose);
    fprintf(outputFile, "vepdef 20 %ld %ld\n", block_size, block_size);
}

void initARC(double centerX, double centerY, double nwaUnit, FILE * outputFile){
    fprintf(outputFile, "[HEADER]\n");
    fprintf(outputFile, "unit_size %0.2f nm/px\n", nwaUnit*1000);
    fprintf(outputFile, "Center_um %0.4f %0.4f\n", centerX, centerY);
    fprintf(outputFile, "Center_px %0.2f %0.2f\n", centerX/nwaUnit, centerY/nwaUnit);
    fprintf(outputFile, "\n\n");
    fprintf(outputFile, "[STITCH]\n\n%d, %4d\n%d, %4d\n\n[DATA]\n\n", 0, 0, 0, 0);
}



void initGDS(FILE * outputFile, unsigned char * gdsPost, unsigned char* polyPre, unsigned char * polyPost, unsigned char * polyForm, int layerNumber){
    int gdspreamble[102] = { 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
        109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
        243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
        0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
        0, 10, 6, 6, 110, 111, 110, 97, 109, 101};
    
    int gdspostamble[8]     = {0, 4, 7, 0, 0, 4, 4, 0};
    
    // {0, 4, 8, 0, 0, 6, 13, 2, 0, [layer = 1], 0, 6, 14, 2, 0, 0};
    int polypreamble[16]    = {0, 4, 8, 0, 0, 6, 13, 2, 0, layerNumber, 0, 6, 14, 2, 0, 0};
    int polypostamble[4]    = {0, 4, 17, 0};
    int polyBlockFormat[4]  = {0, 44, 16, 3};
    
    unsigned char gdsPre[102];
    for (int k = 0; k < 102; k++)
        gdsPre[k] = (unsigned char) gdspreamble[k];
    for (int k = 0; k < 8; k++)
        gdsPost[k] = (unsigned char) gdspostamble[k];
    for (int k = 0; k < 16; k++)
        polyPre[k] = (unsigned char) polypreamble[k];
    for (int k = 0; k < 4; k++)
        polyPost[k] = (unsigned char) polypostamble[k];
    for (int k = 0; k < 4; k++)
        polyForm[k] = (unsigned char) polyBlockFormat[k];

    fwrite(gdsPre, sizeof(char), 102, outputFile);
}

void renderGDS(FILE * outputFile, unsigned char * gdsPost){
    fwrite(gdsPost, sizeof(char), 8, outputFile);
    fclose(outputFile);
}
void renderARC(FILE * outputFile){
    fprintf(outputFile, "[DATA]\n");
    fclose(outputFile);
}
void renderWRV(FILE * outputFile){
    fclose(outputFile);
}

void notifyJoanie(float progress, int zpID){
    char stringBuffer[100];
    for (int k = 0; k < 100; k++)
        stringBuffer[k] = ' ';
    sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%d&progress=%0.3f\"",
            zpID, progress);
    system(stringBuffer);
}


int main(int argc, char** argv)
{
    clock_t begin = clock();

    
    if (argc == 1) {
        printf("No input parameters!!! \n");
        system("PAUSE");
        return 0;
    }
    char** argv_test;
    argv_test = argv + 1;
    double zTol = atof(*(argv_test++));
    double lambda_nm = atof(*(argv_test++));
    double p = atof(*(argv_test++));
    double q = atof(*(argv_test++));

    // Need p and q to be positive
    
    // Need p to define NA, assume NA is always defined on fast side
    int virtualObject = 0;
    if (abs(p) > abs(q)){ // swap
        double tempP = p;
        p = q;
        q = tempP;
    }
    if (p * q < 0){
        virtualObject = -1;
        p = abs(p);
        q = abs(q);

        printf("Detected Virtual swapped");
    }
    
    double obscurationSigma = atof(*(argv_test++));
    double NA = atof(*(argv_test++)); // Changed from p_na to anamorphicI
    int nZerns = atoi(*(argv_test++));
    
    double orders[2 * nZerns];
    for (int k = 0; k < nZerns; k++) {
        orders[k] = atof(*(argv_test++));
        orders[k + nZerns]  = atof(*(argv_test++));
    }

    int customMaskIdx       = atoi(*(argv_test++)); // Changed alpha to custom mask
    double beta             = atof(*(argv_test++));
    double CRAAz            = atof(*(argv_test++))* M_PI/180; // call in degrees
    double CRA              = atof(*(argv_test++)) * M_PI/180;  // call in degrees
    double anamorphicFac    = atof(*(argv_test++));
    double ZPCPhase         = atof(*(argv_test++));
    double APD              = atof(*(argv_test++));
    double APD_window       = atof(*(argv_test++));
    double ZPCR2            = atof(*(argv_test++));
    double ZPCR1            = atof(*(argv_test++));
    double bias_nm          = atof(*(argv_test++));
    int File_format         = atof(*(argv_test++));
    int Opposite_Tone       = atof(*(argv_test++));
    int FSIdx               = atof(*(argv_test++));
    double buttressGapWidth = atof(*(argv_test++));
    double buttressPeriod   = atof(*(argv_test++));
    int setToCenter         = atoi(*(argv_test++));
    //1 million pixel with a pixel size 500 pm (0.5 nm)
    // Thus it transfer to 500 um block size
    long block_size         = atol(*(argv_test++));
    long block_unit_size_pm = 500;
    int NoP                 = atoi(*(argv_test++));
    int IoP                 = atoi(*(argv_test++));
    //int zpID = atoi(*(argv_test++));
    double blockGrid_pm        = atof(*(argv_test++)); // block grid (WRV)


    int layerNumber         = atoi(*(argv_test++));
    bool curl_on            = false;
    int nwaUnitSelection    = atoi(*(argv_test++));
    char * fileName         = *argv_test;
    
    double lambda, bias_um;
    double NA_P = NA + sin(CRA);
    double D = 2*abs(p)*tan(asin(NA_P));
    
    
    
    double nwaUnit = .004; // NWA unit for ARC file in [um/px]
    
    switch (nwaUnitSelection){
        case 0:
            nwaUnit = 1; // not used
            break;
        case 1:
             nwaUnit = 0.00175;
            break;
        case 2:
            nwaUnit = 0.002;
            break;
        case 3:
            nwaUnit = 0.0025;
            break;
        case 4:
            nwaUnit = 0.004;
            break;
        case 5:
            nwaUnit = 0.005;
            break;
        case 6:
            nwaUnit = 0.006;
            break;
        case 7:
            nwaUnit = 0.007;
            break;
        case 8:
            nwaUnit = 0.008;
            break;

            
    }
    
    lambda = lambda_nm / 1000;
    bias_um = bias_nm/1000;

    double blockGrid_pixels = 0;
    // WRV compute blockGrid in units of pixels:
    if (blockGrid_pm >= block_unit_size_pm && blockGrid_pm > 0){
        blockGrid_pixels = blockGrid_pm /   ((double) block_unit_size_pm);
        printf("Rounding WRV pixels by \t\t= %d pixels\n", (int)roundf(blockGrid_pixels));
    } else {
        blockGrid_pixels = 0;
    }
    
    
    printf("\nCreating zone plate:\nNA_P \t\t\t= %0.3f\n", NA_P);
    printf("NA_child \t\t= %f \n", NA);
    printf("lambda \t\t\t= %0.3f nm\n", lambda_nm);
    printf("objective length \t= %g um\n", p);
    printf("image length \t\t= %g um\n", q);
    
    //Aberration order
    for (int j = 0; j < nZerns; j++) {
        printf("Zernike Number \t\t= %f \n", orders[j]);
        printf("Weight \t\t\t= %f waves\n", orders[j + nZerns]);
    }


    printf("WRV block grid \t = %g\n", blockGrid_pm);
    printf("BlockGrid pixels \t = %d\n", (int)roundf(blockGrid_pixels));

    
    printf("tilted angle beta \t= %g in radian\n", beta);
    printf("Max zone phase error\t= lambda/%d\n", (int)(1 / zTol));
    printf("Zone plate diameter \t= %0.2f um\n", D);
    printf("ZPC phase shift \t= %f in radians\n", ZPCPhase);
    printf("Apodized term \t\t= %f %% transmission\n", APD * 100);
    if (APD_window != 0) {
        printf("Apodized window \t= %f \n", APD_window);
    }
    printf("Zone width bias \t\t= %f nm\n", bias_nm);
    printf("Is Free Standing? \t\t= %d \n", FSIdx);
    printf("W\t\t\t= %f dr \n", buttressGapWidth);
    printf("T\t\t\t= %f dr \n", buttressPeriod);
    printf("Block size for wrv file: %ld\n", block_size);
    printf("Generating partial %d of total %d total.\n", IoP, NoP);
    if (Opposite_Tone == 1) {
        printf("Tone swapped");
    }
    int wrv_split;
    if (D >= 500 && File_format == 3 && layerNumber == 1) {
        wrv_split = 1;
        printf("WARNING: WRV is overrunning blocksize");
    }
    printf("\n");
    
    if (File_format == 0) {
        printf("Save as Arc information\n");
    }
    else if (File_format == 1) {
        printf("Save as trapezoid information in txt file\n");
    }
    
    else if (File_format == 2) {
        printf("Save as GDSII file \n");
    }
    else if (File_format == 3) {
        printf("Save as wrv file \n");
    }
    else {
        printf("Save as zone number & dosage information in txt file \n");
    }
    
    if (NoP == 1 && IoP != 1) {
        printf("\nIndex of part has to be 1 if number of parts is 1 !\n");
        system("PAUSE");
        return 0;
        
    }
    if (nwaUnitSelection > 0){
        printf("NWA unit selection is %d with value %0.5f.\n", nwaUnitSelection, nwaUnit);
    }

    
    /*************************/
    // BEGIN NEW ZP CODE:
    
    FILE * outputFile = NULL;
    FILE * supportTextFile = NULL;
    double dbscale = 0; // db unit to microns
    int numBlocksOnSide = 1;
    double rGuess, rGuessp1, Rn, Rnp1, dr, buttressWidth, alphaBT, alphaZT, alpha, 
                x, y, cx, cy, R1, R2, dR1, startAngle, currentAngle, arcStart, phase, 
                RCM = 0, tR1, tR2, f, pNA, RN, rNp, rNq, RNp1, dRN, 
                Rnpa2, Rnma2, tR1pa, tR1ma, tR2pa, tR2ma,
                minDose, maxDose, doseBias, 
                zpCenterX, zpCenterY, offsetX = 0, offsetY = 0, drawAngle;
    double pupilCoordinates[2];

    long totalPoly = 0;
    unsigned char gdsPost[8];
    unsigned char polyPre[16];
    unsigned char polyPost[4];
    unsigned char polyForm[4];
    
    /* Figure out number of zones by counting the number of waves that separate rN and f */
    f           = 1/(1/p + 1/q);
    pNA         = NA + sin(CRA);
    RN          = (p) * tan(asin(pNA));    // R in plane of ZP
    rNp         = sqrt(RN*RN + p*p);     // r is hypotenuse of RN and p
    rNq         = sqrt(RN*RN + q*q);     // r is hypotenuse of RN and p
    zpCenterX   = -(p)*tan(CRA)*sin(CRAAz);
    zpCenterY   = (p)*tan(CRA)*cos(CRAAz);
    
    if (setToCenter){
        offsetX = -zpCenterX;
        offsetY = -zpCenterY;
    }
    if (File_format == 3) { // WRV: shift by one half block size
        numBlocksOnSide = layerNumber; // For WRV files we send num blocks on the GDSLayer input

        offsetX += numBlocksOnSide*block_size*block_unit_size_pm / 2 / 1000000;
        offsetY += numBlocksOnSide*block_size*block_unit_size_pm / 2 / 1000000;
    }

    
    // Compute number of zones in this zone plate using RN path length difference
    int N = 2*(int)((rNp - p + virtualObject*(rNq - q))/lambda);
    printf("Zoneplate has %d zones \n", N);


    // Compute dRN and dR1 to bound dose information
    RN          = secantSolve(sqrt(N*lambda*f), 0, N, p, q, 0, lambda, beta, virtualObject);
    RNp1        = secantSolve(sqrt((N+1)*lambda*f), 0, N + 1, p, q, 0, lambda, beta, virtualObject);
    dRN         = RNp1 - RN;
    R1          = secantSolve(sqrt(lambda*f), 0, 1, p, q, 0, lambda, beta, virtualObject);
    R2          = secantSolve(sqrt(2*lambda*f), 0, 2, p, q, 0, lambda, beta, virtualObject);
    dR1         = R2 - R1;
    maxDose     = dRN/(dRN - bias_um);
    minDose     = dR1/(dR1 - bias_um);
    
    switch (File_format){
        case 0: // NWA (arc)
            if ((outputFile = fopen(strcat(fileName, ".nwa"), "wb")) == NULL)
                printf("Cannot open file.\n");
            initARC(offsetX, offsetY, nwaUnit, outputFile);
            break;
        case 1: // GDS
            dbscale = 10000; // db unit to microns
            if ((outputFile = fopen(strcat(fileName, ".gds"), "wb")) == NULL)
                printf("Cannot open file.\n");
            initGDS(outputFile, gdsPost, polyPre, polyPost, polyForm, layerNumber);
            break;
        case 2://GDS + txt
            dbscale = 10000; // db unit to microns
            if ((outputFile = fopen(strcat(fileName, ".gds"), "wb")) == NULL)
                printf("Cannot open file.\n");
            initGDS(outputFile, gdsPost, polyPre, polyPost, polyForm, layerNumber);
            if ((supportTextFile = fopen(strcat(fileName, ".txt"), "wb")) == NULL)
                printf("Cannot open file.\n");
            
            writeSupportTxtHeader(dR1, dRN, bias_um, supportTextFile);
            break;
        case 3: // WRV
            dbscale = 1000000/block_unit_size_pm; // db unit to microns
            if ((outputFile = fopen(strcat(fileName, ".wrv"), "wb")) == NULL)
                printf("Cannot open file.\n");
            
            initWRV(outputFile, minDose, maxDose, block_size, block_unit_size_pm);
            break;
    }

    long clockSpeed;
    bool isGapZone = false;

    // If we're printing obscuration only, then set N to 0:
    if (customMaskIdx == 16){
        N = 0;
    }

    for (int n = 1; n <= N; n++){
        //Compute initial R at an arbitrary angle.  This will seed information regarding RoC for zone tol:
        rGuess     = sqrt(n*lambda*f);
        rGuessp1   = sqrt((n + 1)*lambda*f);
        Rn         = secantSolve(rGuess, 0, n, p, q, 0, lambda, beta, virtualObject);
        Rnp1       = secantSolve(rGuessp1, 0, n + 1, p, q, 0, lambda, beta, virtualObject);
        dr         = Rnp1 - Rn;
        
        // Write to support text file, may deprecate this shortly
        if(File_format == 2) {
            clockSpeed = getPolyClockSpeed(dR1, dRN, dr, bias_um);
            writeSupportTxtZoneDose(clockSpeed, n, Rn, supportTextFile);
        }
    
        // If buttresses are specified, condition segment width on zone parity
        buttressWidth = 0;
        // ButtressWidth is the width of the buttressed zone segment. Gap is the space in between - technically, the actual buttress
        switch (FSIdx){
            case 0:// % no buttresses
                if (n % 2 == Opposite_Tone){
                    continue;
                }
                buttressWidth = 0;
                break;
            case 1:// % gapped zones
                if (n % 2 == Opposite_Tone){
                    continue;
                }
                buttressWidth = dr*(buttressPeriod - buttressGapWidth);
                break;
            case 2: // full zones + gaps
                if (n % 2 == Opposite_Tone){ //% Even zones get gaps
                    if (n == N) // Don't make last zone of gaps
                        continue;
                    
                    isGapZone = true;
                    buttressWidth = dr*(buttressGapWidth);
                }else{ //% Odd zones get Full zones
                    isGapZone = false;
                    buttressWidth = 0;
                }
                break;
        }
    
        /* We can comply with zone tol constraint by leveraging the radius of curvature of the present zone programmatically. */
        alphaZT = 2*acos(1/(zTol * lambda / Rn + 1)); // Full angle subtended by valid arc segment
        
        // Compare with full angle specified by buttress width
        alphaBT = buttressWidth/Rn;
        
        // Set alpha to the tighter of the two constraints
        if ( (buttressWidth != 0 && alphaBT < alphaZT) ){
            alpha = alphaBT;
        }
        else {
            if (buttressWidth == 0){
                alpha = alphaZT;
            } else {

                // round this to nearest mulltiple
                int alphaRatio = (alphaBT / alphaZT) + 1;
                alpha = alphaBT * 1.01 / ((double) alphaRatio);
            }
        }
        
        // printf("alphaBT: %0.3f, alphaZT: %0.3f, alpha: %0.3f, isGapZone: %d\n", alphaBT, alphaZT, alpha, (int) isGapZone);
        //For dealing with mixed conditions above
        arcStart = -1;
        
        // Log trap count:
        long trapCount = 0;
        
        //Loop through angle
        startAngle     = ((double) (rand() % 1000))/1000 * 2 * M_PI ; // randomize start
        currentAngle = startAngle;

        double gapZoneSize = 0;
        
        while(true){
            if (currentAngle >= startAngle + 2 * M_PI){
                break;
            }
            // Get relative coordinates
            x = Rn*cos(currentAngle);
            y = Rn*sin(currentAngle);
        
            // Convert coordinates to pupil coordinates
            /**
             *  For tilted zone plates we do the transfromation:
             *      x => x cos(beta)
             *      p => p - x sin(beta)
             */
            getPupilCoordinatesFromZPCoordinates(x, y, p, pupilCoordinates, beta, NA, CRA);
            cx = pupilCoordinates[0];
            cy = pupilCoordinates[1];
  
            
            // Accept or reject trap based on pupil boundaries and obscuration
            if (anamorphicFac != 1 && anamorphicFac > 0){ // Anamorphic case
                if (!bIsInAnamorphicPupil(cx, cy, anamorphicFac, obscurationSigma)){
                    currentAngle = currentAngle + alpha;
                    continue;
                }
            } else { // isomorphic case
                if (!bIsInGeometry(cx, cy, obscurationSigma, CRA)){
                    currentAngle = currentAngle + alpha;
                    continue;
                }
            }
            
            
            // Apply custom mask
            if (customMaskIdx != 0){
                if(!bIsInCustomMask(cx, cy, customMaskIdx)){
                    currentAngle = currentAngle + alpha;
                    continue;
                }
            }
        
            // Compute phase terms in waves
            phase = getPhaseTerm(cx, cy, orders, nZerns, ZPCPhase, ZPCR1, ZPCR2);
            // printf("phase: %0.3f,  ", phase);
            phase += customPhase(cx, cy, customMaskIdx)/(2*M_PI);
            // printf("phase: %0.3f,  ", phase);

            // Use initial R as seeds for each subsequent compuptation of zone radii
            Rn   = secantSolve(Rn, currentAngle, n, p, q, phase, lambda, beta, virtualObject);
            Rnp1 = secantSolve(Rnp1, currentAngle, n + 1, p, q, phase, lambda, beta, virtualObject);
            
            // Compute CM to optimized trap to arc and trap coords
            dr  = Rnp1 - Rn;
            RCM = (Rnp1*Rnp1*Rnp1 - Rn*Rn*Rn) / (Rnp1*Rnp1 - Rn*Rn)  * 2/3 * sin(alpha)/alpha;
            //RCM = (Rn + dr/2)*sin(alpha)/alpha; // CM of arc: center trap on arc CM rather than matching
            
            // Solve at +/- alpha/2
            Rnpa2 = secantSolve(Rn, currentAngle + alpha/2, n, p, q, phase, lambda, beta, virtualObject);
            Rnma2 = secantSolve(Rn, currentAngle - alpha/2, n, p, q, phase, lambda, beta, virtualObject);


            // Rotate zone plate by CRA azimuth
            drawAngle = currentAngle + CRAAz;
            
            if (File_format == 0){ // ARC format
                doseBias = dr/(dr - bias_um); // Inverse of zone area bias
                if (isGapZone){
                    exportArc(Rn - bias_um/2, dr + bias_um, drawAngle, alpha, doseBias, nwaUnit, offsetX, offsetY, outputFile);
                } else {
                    exportArc(Rn + bias_um/2, dr - bias_um, drawAngle, alpha, doseBias, nwaUnit, offsetX, offsetY, outputFile);
                }
            }
            else { // GDS or WRV format
                
                
                
                if (isGapZone){
                    tR1 = (RCM - dr/2)*1/cos(alpha/2) - bias_um/2;
                    tR2 = (RCM + dr/2)*1/cos(alpha/2) + bias_um/2;
                } else {
                    tR1 = (RCM - dr/2)*1/cos(alpha/2) + bias_um/2;
                    tR2 = (RCM + dr/2)*1/cos(alpha/2) - bias_um/2;
                }

                // Compute corrected radii at the edges of trap instead of center:
                tR1ma = tR1 + (Rnma2 - Rn);
                tR1pa = tR1 + (Rnpa2 - Rn);
                tR2ma = tR2 + (Rnma2 - Rn);
                tR2pa = tR2 + (Rnpa2 - Rn);
  
                // Coordinates of trap
                double trapCoords_um[8];
                    trapCoords_um[0] = (tR1ma*cos(drawAngle - alpha/2) + offsetX);
                    trapCoords_um[1] = (tR1ma*sin(drawAngle - alpha/2) + offsetY);
                    trapCoords_um[2] = (tR1pa*cos(drawAngle + alpha/2) + offsetX);
                    trapCoords_um[3] = (tR1pa*sin(drawAngle + alpha/2) + offsetY);
                    trapCoords_um[4] = (tR2pa*cos(drawAngle + alpha/2) + offsetX);
                    trapCoords_um[5] = (tR2pa*sin(drawAngle + alpha/2) + offsetY);
                    trapCoords_um[6] = (tR2ma*cos(drawAngle - alpha/2) + offsetX);
                    trapCoords_um[7] = (tR2ma*sin(drawAngle - alpha/2) + offsetY);
            
                // Coordinates rounded to pixel scale
                long trapCoords[10];
                    trapCoords[0] = (long) dbscale*trapCoords_um[0];
                    trapCoords[1] = (long) dbscale*trapCoords_um[1];
                    trapCoords[2] = (long) dbscale*trapCoords_um[2];
                    trapCoords[3] = (long) dbscale*trapCoords_um[3];
                    trapCoords[4] = (long) dbscale*trapCoords_um[4];
                    trapCoords[5] = (long) dbscale*trapCoords_um[5];
                    trapCoords[6] = (long) dbscale*trapCoords_um[6];
                    trapCoords[7] = (long) dbscale*trapCoords_um[7];
                    trapCoords[8] = (long) dbscale*trapCoords_um[0];
                    trapCoords[9] = (long) dbscale*trapCoords_um[1];


                // compute clockspeed by area fraction:
                // clockSpeed = getPolyClockSpeedFromAreaFraction(1/maxDose, 1/minDose, trapCoords_um, alpha * RCM * dr, dbscale, blockGrid_pixels);
                clockSpeed = getPolyClockSpeed(dR1, dRN, dr, bias_um);
            
                // Export shape
                exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format, clockSpeed, blockGrid_pixels);
            }
        
            trapCount++;
        
            // Increment angle by amount depending on whether we require a gap or another trap to satisfy zTol requirement
            if (buttressWidth == 0){
                currentAngle = currentAngle + alpha;
            } else if (isGapZone){
                // If the bridges in the odd zones require more than one trap to make, then keep a running tab of size
                gapZoneSize += alpha;
                if (gapZoneSize >= alphaBT){
                    currentAngle = currentAngle + dr*buttressPeriod/RCM - alphaBT;
                    gapZoneSize = 0;
                } else {
                    currentAngle = currentAngle + alpha;
                }
                
            } else{
                if (alphaBT > alphaZT){
                    if (arcStart < 0){ // This segment requires multipe traps to satsify ztol.  Start counting
                        arcStart = currentAngle; 
                    }
                    if ((currentAngle - arcStart + alpha) > alphaBT){ // Require a gap here, reset counter
                        currentAngle = currentAngle + alpha + dr*buttressGapWidth/RCM;
                        arcStart = -1;
                    }
                    else{ // keep adding traps in this segment
                        currentAngle = currentAngle + alpha;
                    }
                }
                else{
                    currentAngle = currentAngle + dr*buttressPeriod/RCM; // This segment satisfies ztol
                }
            }
            
        } // End angle loop
            
        totalPoly += trapCount;
        if (File_format == 0){ // ARC format
            printf("Finished zone %d with %ld arcs\n", n, trapCount);
        } else {
            printf("Finished zone %d with %ld traps.  \tR_%d = %0.3f \n", n, trapCount, n, Rn);
        }
        if (curl_on && n % 100 == 0) {
            float progress = ((float)n)/((float)N);
            notifyJoanie(progress, 1);
        }
        
    } // End zone loop


    // Write obscuration.  Note: cannot write obscuration in arcs presently
    if (Opposite_Tone && obscurationSigma > 0 && File_format != 0){
        printf("Writing obscuration with sigma: %0.2f\n", obscurationSigma);
        printf("Obscuration anamorphic factor: %0.2f\n", anamorphicFac);
        if (anamorphicFac != 1){
            printf("Writing asymmetric obscuration\n");
        }


        double kx, ky;
        int nObsPts = 360;
        double qXCoords[nObsPts];
        double qYCoords[nObsPts];

        double coordinates[2];

        /**
         * Refer to overleaf notes "Generalized zone plate aperture geometry"
         */

        double thStep = 2*M_PI/nObsPts;
        double theta;
        for (int k = 0; k < nObsPts; k++){
            theta = thStep * k + M_PI/2;

            kx = NA*obscurationSigma*cos(theta);
            ky = NA*obscurationSigma*sin(theta)/anamorphicFac + sin(CRA);

            getZPCoordinatesFromFrequency(kx, ky, p, coordinates, beta);

            qXCoords[k] = coordinates[0];
            qYCoords[k] = coordinates[1];
        }

        double azRot = CRAAz;
        for (int k = 0; k < (nObsPts/2 - 1); k++){
            totalPoly++;

            long trapCoords[10];
            trapCoords[0] = (long) dbscale*(offsetX + cos(azRot)*qXCoords[k] - sin(azRot) * qYCoords[k]);
            trapCoords[1] = (long) dbscale*(offsetY + sin(azRot)*qXCoords[k] + cos(azRot) * qYCoords[k]);

            trapCoords[2] = (long) dbscale*(offsetX + cos(azRot)*qXCoords[k+1] - sin(azRot) * qYCoords[k+1]);
            trapCoords[3] = (long) dbscale*(offsetY + sin(azRot)*qXCoords[k+1] + cos(azRot) * qYCoords[k+1]);

            trapCoords[4] = (long) dbscale*(offsetX + cos(azRot)*qXCoords[nObsPts - k - 2] - sin(azRot) * qYCoords[nObsPts - k - 2]);
            trapCoords[5] = (long) dbscale*(offsetY + sin(azRot)*qXCoords[nObsPts - k - 2] + cos(azRot) * qYCoords[nObsPts - k - 2]);

            trapCoords[6] = (long) dbscale*(offsetX + cos(azRot)*qXCoords[nObsPts - k - 1] - sin(azRot) * qYCoords[nObsPts - k - 1]);
            trapCoords[7] = (long) dbscale*(offsetY + sin(azRot)*qXCoords[nObsPts - k - 1] + cos(azRot) * qYCoords[nObsPts - k - 1]);

            trapCoords[8] = (long) dbscale*(offsetX + cos(azRot)*qXCoords[k] - sin(azRot) * qYCoords[k]);
            trapCoords[9] = (long) dbscale*(offsetY + sin(azRot)*qXCoords[k] + cos(azRot) * qYCoords[k]);
        
            // Export shape
            exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format, 65535, blockGrid_pixels);

        }
    }
 
    
    switch (File_format){
        case 0: // arc
            renderARC(outputFile);
            break;
        case 1: // GDS
            renderGDS(outputFile, gdsPost);
            break;
        case 2: // GDS + txt
            renderGDS(outputFile, gdsPost);
            fclose(supportTextFile);
            break;
        case 3: // WRV
            renderWRV(outputFile);
            break;
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Finished zone plate with %ld shapes\n", totalPoly);
    printf("\nZP Generation took: %0.3f seconds\n", elapsed_secs);
    
    if (curl_on)
        notifyJoanie(1, 1);
    
    return 0;
}




