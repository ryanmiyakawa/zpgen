

#include <iostream>
#include <stdio.h>
#include "Eigen/Dense"
#include "Eigen/Core"
#include <fstream>
#include <cmath>
#include <array>
#include <algorithm>
#include <valarray>
#include <string.h>
#include <errno.h>
#include <ctime>
#include <vector>
#include <cstdint>


using namespace std;
using Eigen::MatrixXd;
using Eigen::Matrix3i;
using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::VectorXf;
using Eigen::Dynamic;
using Eigen::ArrayBase;
using Eigen::Array3d;





// Recursively computes the binomial coefficient nCk:\

double nChooseK(double n, double k) {
    if (k == 0 && n >= 0)
        return 1;
    else if (n == 0 && k > 0)
        return 0;
    else
        return nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
}


// Retrieves the value the Zernike Polynomial order J at the coordinate [r, th]:\

double zgenpt(int j, double r, double th) {
    
    // Get dual index [n,m] from j:\
    
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
    
    
    // Compute value of radial function:\
    
    int p = (n - abs(m)) / 2;
    double kParity = 1;
    double Rv = 0;
    
    for (int k = 0; k <= p; k += 1) {
        if (k % 2 == 1)
            kParity = -1;
        
        Rv += kParity*nChooseK((double)(n - k), (double)k)
        *nChooseK((double)(n - 2 * k), (double)(p - k))
        *pow(r, (double)(n - 2 * k));
    }
    
    
    // Compute value of azimuthal function:\
    
    double Av = 0;
    if (m > 0)
        Av = cos(m*th);
    else if (m < 0)
        Av = -sin(m*th);
    else
        Av = 1;
    
    // Combine:\
    
    return Rv*Av;
    
}


double objectiveFn(double r, double th, double * orders, int rowLen,
                   double zoneNum, double p, double lambda, double q,
                   double zpRadius, double PC_term, int PC_on, int offaxis_aberration,
                   double CRA, double NA_C, double beta) {
    
    
    
    return 0;
}


// Finds the zero of a function using the secant method, which is a
// discretized version of Newton's method
// Use p & q instead of Mag & focus
// Include tilting (alpha, beta) into consideration
// In this function, we don't have to use xm1 & xm2 to do the iteration,
// we have to calculate the new radius (xm1?) then put it back to iteration.

double secantSolve(double xm1, double theta, double * orders, int rowLen,
                   double zoneNum, double p, double lambda, double q,
                   double alpha, double beta, double gamma, double zpRadius,
                   int currentIter, int maxIter, double tolX, double PC_term,
                   int PC_on, int offaxis_aberration, double CRA, double NA_C) {
    
    
        return 0;
}



double APD_WindowF(double r, double th, double * orders, int rowLen,
                   double zoneNum, double p, double lambda, double q,
                   double zpRadius, double PC_term, int PC_on, int offaxis_aberration,
                   double CRA, double NA_C, double APD_window) {
    //This function is used to check the APD value at each position relative to the center of the off-axis zoneplate
    //APD_window == 1: Hamming window
    //APD_window == 2: Gaussian window
    
    double APD_value;
    
    double x_0 = r*cos(th);
    double y_0 = r*sin(th);
    
    double Rc = p*tan(CRA / 180 * M_PI);
    double y_offaxis = y_0 - Rc;
    
    double phi = atan2(y_offaxis, x_0);
    th = phi;
    
    // Ken's equation has his starting point at y-axis when phi = 0
    phi = phi - M_PI / 2;
    
    double x_boundary = (p*(NA_C*sin(phi) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(phi) - NA_C*NA_C)));
    double y_boundary = (p*((sin((CRA / 180)*M_PI) + NA_C*cos(phi)) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(phi) - NA_C*NA_C)));
    
    double r_new = sqrt(x_0*x_0 + y_offaxis*y_offaxis);
    zpRadius = sqrt(x_boundary*x_boundary + (y_boundary - Rc)* (y_boundary - Rc));
    
    double rn = r_new / zpRadius;
    //double sigma = 2;
    
    //With rn and th, you can find the APD_value based on it
    //Wait for the equation
    
    if (APD_window == 1) {
        APD_value = 0.54 + (1 - 0.54) * cos(M_PI * rn);
    }
    else if (APD_window == 2) {
        APD_value = exp(-rn*rn / (2 * 1 * 1));
    }
    else if (APD_window == 3) {
        APD_value = exp(-rn*rn / (2 * 2 * 2));
    }
    
    return APD_value;
    
}
// This is a function for CAD_tool only, to determine the block number for wrv file
int BlockNumberFinder(float* xArray, int block_size)
{
    // This is the boundary of the central block
    // and we only use 80% of it: 400 um out of 500 um.
    int new_center = block_size / 2 / 2;
    int boundary = (block_size / 2)*0.8;
    
    VectorXf vertices_coord(8);
    VectorXf vertices_coord3(8);
    
    for (int k = 0; k < 1; k++) {
        
        //printf("%lld out of total trap %lld \n", k+1, numTrap+1);
        //printf("Total polygon %ld \n", totalPoly);
        
        //Turn the coordinates from um to A.
        for (int m = 0; m < 8; m++) {
            vertices_coord(m) = xArray[1 * m + k] * 10000;
        }
        
        //Do the rotation first: -90 degree
        vertices_coord3(0) = vertices_coord(1) + new_center;
        vertices_coord3(1) = -1 * vertices_coord(0) + new_center;
        vertices_coord3(2) = vertices_coord(3) + new_center;
        vertices_coord3(3) = -1 * vertices_coord(2) + new_center;
        vertices_coord3(4) = vertices_coord(5) + new_center;
        vertices_coord3(5) = -1 * vertices_coord(4) + new_center;
        vertices_coord3(6) = vertices_coord(7) + new_center;
        vertices_coord3(7) = -1 * vertices_coord(6) + new_center;
    }
    //The following section is used to determined the belonging of the trapezoide.
    // Using the boundary of the center block as the checking point:
    // x < 50 & y > 450 -> zone 1
    // x > 450 & y < 50 -> zone 9
    VectorXd x_idx(4);
    VectorXd y_idx(4);
    int isFive = 0;
    
    Matrix3i block_number_m;
    block_number_m << 1, 2, 3,
    4, 5, 6,
    7, 8, 9;
    VectorXd trap_blockN(4);
    int most = 0;
    int toReturn = 0;
    
    //Sort out the x & y position
    for (int idx = 0; idx < 8; idx = idx + 2) {
        float test_x = vertices_coord3(idx);
        if (test_x > 4500000) {
            x_idx(idx / 2) = 2;
        }
        else if (test_x < 500000) {
            x_idx(idx / 2) = 0;
        }
        else {
            x_idx(idx / 2) = 1;
        }
    }
    for (int idx = 1; idx < 8; idx = idx + 2) {
        float test_y = vertices_coord3(idx);
        if (test_y > 4500000) {
            y_idx(floor(idx / 2)) = 0;
        }
        else if (test_y < 500000) {
            y_idx(floor(idx / 2)) = 2;
        }
        else {
            y_idx(floor(idx / 2)) = 1;
        }
    }
    // Find out the block number of each point of the trapezoid
    for (int idx = 0; idx < 4; idx = idx + 1) {
        trap_blockN(idx) = block_number_m(y_idx(idx), x_idx(idx));
        //Only need to spot 1 point in zone 5
        if (trap_blockN(idx) == 5) {
            isFive = 1;
        }
    }
    
    //Sort out the block number distribution
    VectorXi blockN_count = VectorXi::Zero(9);
    
    //1st element starts at 0!
    
    for (int idx = 0; idx < 4; idx = idx + 1) {
        blockN_count(trap_blockN(idx) - 1) = blockN_count(trap_blockN(idx) - 1) + 1;
    }
    
    //Block number rule:
    //Every point is in the same zone -> use the zone number
    //Any point is in block 5 -> block_number =5
    //choose the 1 has most point (at least 2) in it
    //if 2 vs. 2 choose	smaller zone nunmber
    std::ptrdiff_t max_idx;
    
    if (blockN_count.maxCoeff(&max_idx) >= 4) {
        return (max_idx + 1);
    }
    else if (isFive == 1) {
        return 5;
    }
    else {
        //Find out which zone number has most points in it
        for (int idx = 0; idx < 9; idx = idx + 1) {
            if (blockN_count(idx) > most) {
                most = blockN_count(idx);
                toReturn = idx;
            }
        }
        return toReturn + 1;
    }
}

int nfinder(double R0, double alpha, double beta, double gamma, double lambda, double p, double q)
{
    VectorXd res(1001);
    int i = 0;
    double x_max;
    double y_max;
    double x_out;
    double y_out;
    double z_out;
    
    // We have to include tilt angle here with coordinate transform to get the right ns & ne
    for (double theta = 0; theta <= 2 * M_PI; theta = (2 * M_PI) / 1000 + theta)
    {
        x_max = R0* cos(theta);
        y_max = R0* sin(theta);
        
        //Instead of doing matrix multiplication, we get the result of each element.
        
        //x_out = x_max*cos(beta);
        //y_out = x_max*1*sin(alpha)*sin(beta) + y_max*cos(alpha);
        //z_out = x_max*-1*cos(alpha)*sin(beta) + y_max*sin(alpha);
        
        //x_out = x_max* cos(beta)* cos(gamma) - y_max*cos(beta)*sin(gamma);
        //y_out = x_max*(sin(alpha)* sin(beta)* cos(gamma)+ cos(alpha)* sin(gamma))+ y_max* (cos(alpha)* cos(gamma)- sin(alpha)*sin(beta)*sin(gamma));
        //z_out = x_max*(sin(alpha)* sin(gamma)-cos(alpha)*sin(beta)*cos(gamma))+ y_max*(sin(alpha)*cos(gamma)+ cos(alpha)*sin(beta)*sin(gamma));
        
        
        // find out the angle between vector and +z axis (001)
        
        //double numerator = z_out*1;
        //double denominator =  sqrt(x_out*x_out + y_out*y_out +z_out*z_out)*1;
        
        //double th = acos( numerator/ denominator);
        double th = 1.5708;
        
        //calculate OPD
        
        res(i) = sqrt(p*p + R0*R0 + 2 * p*R0*cos(th)) - p + sqrt(q*q + R0*R0 - 2 * q*R0*cos(th)) - q;
        i++;
        
    }
    
    int n = (int)(res.minCoeff() * 2 / lambda);
    return n;
}

double* CoordinateTransform(double R2, double alpha, double beta, double gamma,
                            int numpoints, int & CT_array_size, double PC_term,
                            int OuterBoundary, double NA_C, double p, double q
                            , double CRA, double APD) {
    // double R1 was used to define the center of the offaxis zone plate
    // Now we use p*tan(CRA) to define Rc
    
    // double Rc = (p*tan(CRA);
    
    // Define points on ZP plane
    
    MatrixXd xyz_plane(3, numpoints);
    
    VectorXd th(numpoints);
    
    double f = (p*q) / (p + q);
    
    double Rc = p*tan((CRA / 180)*M_PI);
    
    double x_test, y_test, R_test;
    
    for (int n = 0; n < numpoints; n++) {
        th(n) = (2 * M_PI / numpoints)*n;
    }
    if ((PC_term != 0 || APD != 1) && OuterBoundary == 0) {
        
        for (int i = 0; i < numpoints; i++) {
            // This part is used to define circular shape
            // theta is not needed anymore, use gamma to rotate it later.
            // Need to find a way to define the center of the zone plate by Ken's equation
            // and move the phase shifted regionto match the zone plate region
            // Center point (0, Rc)
            
            xyz_plane(0, i) = ((R2*R2) / sqrt((R2*cos(th(i))*R2*cos(th(i))) + (R2*sin(th(i))*R2*sin(th(i))))) * cos(th(i));
            xyz_plane(1, i) = ((R2*R2) / sqrt((R2*cos(th(i))*R2*cos(th(i))) + (R2*sin(th(i))*R2*sin(th(i))))) * sin(th(i)) + Rc;
            xyz_plane(2, i) = 0;
            
        }
    }
    else {
        for (int i = 0; i < numpoints; i++) {
            
            //xyz_plane(0, i) = ((R2*R2*2)/sqrt((R2*2*cos(th(i))*R2*2*cos(th(i))) + (R2*sin(th(i))*R2*sin(th(i))))) * cos(th(i)) + R1*cos(theta);
            //xyz_plane(1, i) = ((R2*R2*2)/sqrt((R2*2*cos(th(i))*R2*2*cos(th(i))) + (R2*sin(th(i))*R2*sin(th(i))))) * sin(th(i)) + R1*sin(theta);
            //xyz_plane(2, i)	= 0;
            xyz_plane(0, i) = (p*(NA_C*sin(th(i)) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(th(i)) - NA_C*NA_C)));
            xyz_plane(1, i) = (p*((sin((CRA / 180)*M_PI) + NA_C*cos(th(i))) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(th(i)) - NA_C*NA_C)));
            xyz_plane(2, i) = 0;
            
            
        }
    }
    
    
    // Tilt it follow yz xz order
    // Alpha: rotate angle alpha about x axis
    // Beta: rotate angle beta about y axis
    
    MatrixXd uvw_plane(3, numpoints);
    
    MatrixXd alphaTransform(3, 3);
    
    alphaTransform(0, 0) = 1;
    alphaTransform(0, 1) = 0;
    alphaTransform(0, 2) = 0;
    
    alphaTransform(1, 0) = 0;
    alphaTransform(1, 1) = cos(alpha);
    alphaTransform(1, 2) = -sin(alpha);
    
    alphaTransform(2, 0) = 0;
    alphaTransform(2, 1) = sin(alpha);
    alphaTransform(2, 2) = cos(alpha);
    
    MatrixXd betaTransform(3, 3);
    
    betaTransform(0, 0) = cos(beta);
    betaTransform(0, 1) = 0;
    betaTransform(0, 2) = sin(beta);
    
    betaTransform(1, 0) = 0;
    betaTransform(1, 1) = 1;
    betaTransform(1, 2) = 0;
    
    betaTransform(2, 0) = -sin(beta);
    betaTransform(2, 1) = 0;
    betaTransform(2, 2) = cos(beta);
    
    
    MatrixXd gammaTransform(3, 3);
    
    gammaTransform(0, 0) = cos(gamma);
    gammaTransform(0, 1) = -sin(gamma);
    gammaTransform(0, 2) = 0;
    
    gammaTransform(1, 0) = sin(gamma);
    gammaTransform(1, 1) = cos(gamma);
    gammaTransform(1, 2) = 0;
    
    gammaTransform(2, 0) = 0;
    gammaTransform(2, 1) = 0;
    gammaTransform(2, 2) = 1;
    
    
    // Matrix multiplication to finish the transform
    
    uvw_plane = alphaTransform* betaTransform* gammaTransform* xyz_plane;
    
    //Transfer it from uvw_plane to CT_array
    
    
    VectorXd CT_array0(numpoints * 3);
    
    for (int j = 0; j < numpoints; j++) {
        
        CT_array0(0 + j * 3) = uvw_plane(0, j);
        CT_array0(1 + j * 3) = uvw_plane(1, j);
        CT_array0(2 + j * 3) = uvw_plane(2, j);
    }
    
    CT_array_size = numpoints * 3;
    double* CT_array = new double[CT_array_size];
    
    
    for (int i = 0; i < CT_array_size; i++) {
        
        CT_array[i] = CT_array0(i);
        
    }
    
    return CT_array;
}

double* define_function(double* CT_array, double p, double q, double lambda, double NA_P, double alpha,
                        double beta, double gamma, double zpRadius, double * orders, int rowLen, int CT_array_size,
                        int & nandtheta_array_size, int &ZoneNum_min, int &ZoneNum_max, double offaxis_aberration, double CRA, double NA_C) {
    
    int size = CT_array_size / 3;
    MatrixXd nandtheta(2, size);
    
    for (int i = 0; i < size; i++) {
        
        double x = CT_array[0 + i * 3];
        double y = CT_array[1 + i * 3];
        
        double R = sqrt(x*x + y*y);
        double th;
        double ZoneNum;
        
        if (x >= 0 && y >= 0) {
            th = asin(y / R);
        }
        else if (x < 0 && y >= 0) {
            th = M_PI - asin(y / R);
        }
        else if (x <= 0 && y<0) {
            th = M_PI + asin(abs(y / R));
        }
        else {
            th = 2 * M_PI - asin(abs(y / R));
        }
        
        ZoneNum = secantSolve2(R, th, orders, rowLen, p, lambda, q, alpha, beta, gamma, zpRadius, offaxis_aberration, CRA, NA_C);
        
        nandtheta(0, i) = ZoneNum;
        nandtheta(1, i) = th;
        
    }// End of for loop
    
    VectorXd ZoneNum = nandtheta.row(0);
    
    double test0 = ZoneNum.minCoeff();
    double test1 = ZoneNum.maxCoeff();
    
    ZoneNum_min = (int)ceil(ZoneNum.minCoeff());
    ZoneNum_max = (int)floor(ZoneNum.maxCoeff());
    
    printf("off-axis starting zone number = %d\n", ZoneNum_min);
    printf("off-axis ending zone number = %d\n", ZoneNum_max);
    
    //Find the crossing point for two non-integer Zonenumber & and get the angle of it.
    
    MatrixXd nandtheta_test(2, 1);
    int l = 0;
    int MatrixSize = 1;
    
    for (int i = 0; i < (ZoneNum.size()); i++) {
        
        double a;
        double b;
        double th_a;
        double th_b;
        
        if (i == ZoneNum.size() - 1) {
            a = nandtheta(0, i);
            b = nandtheta(0, 0);
            
            th_a = nandtheta(1, i);
            th_b = nandtheta(1, 0);
        }
        else {
            a = nandtheta(0, i);
            b = nandtheta(0, i + 1);
            th_a = nandtheta(1, i);
            th_b = nandtheta(1, i + 1);
        }
        if (a>b) {
            a = (int)floor(a);
            b = (int)ceil(b);
        }
        else {
            a = (int)ceil(a);
            b = (int)floor(b);
        }
        
        if (a == b) {
            
            int n = (int)a;
            nandtheta_test(0, l) = n;
            // Consider the angle. Solve the issue when it cross 2*PI
            if (min(th_a, th_b)<M_PI / 2 && max(th_a, th_b)>3 * M_PI / 2)
            {
                nandtheta_test(1, l) = ((th_a + 2 * M_PI) + th_b) / 2;
                if (nandtheta_test(1, l) > 2 * M_PI)
                    nandtheta_test(1, l) = nandtheta_test(1, l) - 2 * M_PI;
            }
            else {
                nandtheta_test(1, l) = (th_a + th_b) / 2;
            }
            
            l++;
            MatrixSize++;
            nandtheta_test.conservativeResize(2, MatrixSize);
        }//End of a==b if loop;
        
    }//End of for loop
    
    nandtheta_test.conservativeResize(2, --MatrixSize);
    
    
    
    // This matrix is to sort out the ZoneNum, and their starting and ending angle
    //int size_test= nandtheta_test.size;
    MatrixXd nandtheta2(3, ZoneNum_max - ZoneNum_min + 1);
    int index = 0;
    
    for (int k = ZoneNum_min; k < (ZoneNum_max + 1); k++) {
        
        int i = 0;// for nandtheta_test index
        int j = 0;// to check wheather it is 1st time or not
        
        while (1) {
            if (nandtheta_test(0, i) == k) {
                
                if (j == 0) {
                    nandtheta2(0, index) = nandtheta_test(0, i);
                    nandtheta2(1, index) = nandtheta_test(1, i);
                    
                    j++;
                }// End of j = 0 loop
                else {
                    nandtheta2(2, index) = nandtheta_test(1, i);
                    break; // break while loop
                }
            }// End of if loop
            
            i++; // keep searching for the second one
        }// End of while loop
        
        
        //put the min theta at 2nd row
        
        double min_theta = min(nandtheta2(1, index), nandtheta2(2, index));
        double max_theta = max(nandtheta2(1, index), nandtheta2(2, index));
        
        nandtheta2(1, index) = min_theta;
        nandtheta2(2, index) = max_theta;
        
        
        index++;
    }//End of nandtheta2 for loop
    
    
    double* nandtheta_array = new  double[(ZoneNum_max - ZoneNum_min + 1) * 3];
    
    //Save all the information into array
    
    for (int j = 0; j < (ZoneNum_max - ZoneNum_min + 1); j++) {
        
        nandtheta_array[0 + j * 3] = nandtheta2(0, j);
        nandtheta_array[1 + j * 3] = nandtheta2(1, j);
        nandtheta_array[2 + j * 3] = nandtheta2(2, j);
    }
    
    nandtheta_array_size = (ZoneNum_max - ZoneNum_min + 1) * 3;
    return nandtheta_array;
}

//This function is used to check if the point is within the off-axis zonepalte
// output "1" if the point is inside the boundary.
// output "0" if the point is out of the boundary.
// Reuse the code for off-axis aberration to determine if the Radisu from the cneter of the
// off-axis zoneplate is larger than 1.

int IsOffAxis(double r, double th, double p, double q, double CRA, double NA_C)
{
    
    
    double x_0, y_0, Rc, y_offaxis, r_new, zpRadius, rn, phi, x_boundary, y_boundary;
    
    x_0 = r*cos(th);
    y_0 = r*sin(th);
    
    Rc = p*tan(CRA / 180 * M_PI);
    y_offaxis = y_0 - Rc;
    
    phi = atan2(y_offaxis, x_0);
    th = phi;
    
    // Ken's equation has his starting point at y-axis when phi = 0
    phi = phi - M_PI / 2;
    
    x_boundary = (p*(NA_C*sin(phi) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(phi) - NA_C*NA_C)));
    y_boundary = (p*((sin((CRA / 180)*M_PI) + NA_C*cos(phi)) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(phi) - NA_C*NA_C)));
    
    r_new = sqrt(x_0*x_0 + y_offaxis*y_offaxis);
    zpRadius = sqrt(x_boundary*x_boundary + (y_boundary - Rc)* (y_boundary - Rc));
    
    rn = r_new / zpRadius;
    
    if (rn <= 1) {
        return 1;
    }
    else {
        return 0;
    }
}

// This function is used to generate off-axis ZP

float * Control(double* nandtheta_array, double* nandtheta_array2, double* nandtheta_array3, double lambda, double p, double q, double alpha, double beta, double gamma,
                int ns, int ne, double obscurationSigma, double InnerSigma, double OuterSigma, double NA_P, double* orders, int numA, int &totalArcs,
                double zTol, int nandtheta_array_size, int nandtheta_array_size2, int nandtheta_array_size3, double PC_term, double APD, double bias, int32_t zpID,
                int offaxis_aberration, double CRA, double NA_C, int file_format, char** argv_test, double NoP, double IoP, double APD_window, int curl_on) {
    
    
    double NA_Pt = tan(asin(NA_P));
    double D = 2 * (p*q / (p + q))*NA_Pt;
    double R0 = D / 2;
    double r0 = R0* obscurationSigma;
    int Nmin = 0;
    int	Nmin2 = 0;
    int Nmin3 = 0;
    int Nmax = 0;
    int Nmax2 = 0;
    int Nmax3 = 0;
    int delta = 0;
    int delta2 = 0;
    int delta3 = 0;
    int Nmin0 = 0;
    int Nmin02 = 0;
    int Nmin03 = 0;
    int Nmax0 = 0;
    int Nmax02 = 0;
    int Nmax03 = 0;
    int State = 0;
    double Rc = p*tan((CRA / 180.0)*M_PI);
    
    // For scale bar
    char stringBuffer[100];
    
    
    // sort out the situation we have for following caculation
    
    // Only off-axis
    if (PC_term == 0) {
        State = 1;
    }
    //Only PC
    if (PC_term != 0 && APD == 1) {
        // Only Disk PC
        if (InnerSigma == 0) {
            State = 2;
        }
        else {
            // Annular & PC
            State = 3;
        }
    }
    // PC & APD
    if (PC_term != 0 && APD != 1) {
        // PC & APD & Disk
        if (InnerSigma == 0) {
            State = 4;
        }
        else {
            // PC & APD & Annular
            State = 5;
        }
        
        
    }
    //Only APD
    if (PC_term == 0 && APD != 1) {
        if (APD != 0) {
            // Only Disk APD
            if (InnerSigma == 0) {
                State = 6;
            }
            else {
                // Annular & APD
                State = 7;
            }
        }
        else {
            // Dark field disk
            if (InnerSigma == 0) {
                State = 8;
            }
            else {
                // Dark field annular
                State = 9;
            }
        }
    }
    
    if (APD_window != 0) {
        // Specific APD window application
        State = 10;
    }
    
    
    // Compare to the index we find based on the shape we are interested
    // Extract ZoneNum from nandtheta_array
    VectorXd zoneNum(nandtheta_array_size / 3);
    VectorXd theta_S(nandtheta_array_size / 3);
    VectorXd theta_E(nandtheta_array_size / 3);
    
    VectorXd zoneNum2(nandtheta_array_size2 / 3);
    VectorXd theta_S2(nandtheta_array_size2 / 3);
    VectorXd	theta_E2(nandtheta_array_size2 / 3);
    
    VectorXd zoneNum3(nandtheta_array_size3 / 3);
    VectorXd theta_S3(nandtheta_array_size3 / 3);
    VectorXd	theta_E3(nandtheta_array_size3 / 3);
    
    
    for (int i = 0; i < nandtheta_array_size / 3; i++) {
        
        zoneNum(i) = nandtheta_array[0 + i * 3];
        theta_S(i) = nandtheta_array[1 + i * 3];
        theta_E(i) = nandtheta_array[2 + i * 3];
    } // End of for loop
    
    Nmin0 = (int)zoneNum.minCoeff();
    Nmax0 = (int)zoneNum.maxCoeff();
    
    //Define the "real" Nstart & Nend
    if (ns > Nmin0) {
        Nmin = ns;
        delta = ns - Nmin0;
    }
    else {
        Nmin = Nmin0;
        delta = 0;
    }
    
    if (ne < Nmax0) {
        Nmax = ne;
    }
    else {
        Nmax = Nmax0;
    }
    
    // Make sure it start from the same even/odd order as non off-axis case
    
    if ((ns % 2) != 1) {
        if ((Nmin % 2) == 1) {
            Nmin = Nmin + 1;
            delta = delta + 1;
        }
    }//End of if loop
    else {
        if ((Nmin % 2) != 1) {
            Nmin = Nmin + 1;
            delta = delta + 1;
        }
    }// End of else loop
    
    // if Nmax is controled by Nmax0, we have to -1 becasue we don't have the
    // information about the starting and ending angle of that zone.
    
    if (((Nmax - Nmin + 1) % 2) == 1) {
        if (Nmax == Nmax0) {
            Nmax = Nmax - 1;
        }
        else {
            Nmax = Nmax + 1;
        }
    }//End of if loop
    
    if (Nmin >= Nmax) {
        printf("\n Nmin > Nmax, no ZP allowed\n");
        system("PAUSE");
        return 0;
    }
    
    // Reproduce ZP: collect thetaS & thetaE first for each zone.
    // Then call it into HWcoustomZP as an array to reproduce zp
    // Nmin & Nmax is the starting and ending point for each circle
    int m = 0;
    int m2 = 0;
    int m3 = 0;
    int Ns = delta;
    int Ns2 = delta2;
    int Ns3 = delta3;
    int Ne = Ns + 1;
    int Ne2 = Ns2 + 1;
    int Ne3 = Ns3 + 1;
    int Ns_index = Nmin;
    int Ns_index2 = Nmin2;
    int Ns_index3 = Nmin3;
    int situation = 0;
    int situation2 = 0;
    int situation3 = 0;
    double thetaS0, thetaS02, thetaS03, thetaS0_1, thetaS0_2, thetaS0_1_2, thetaS0_2_2, thetaS0_1_3, thetaS0_2_3;
    double thetaE0, thetaE02, thetaE03, thetaE0_1, thetaE0_2, thetaE0_1_2, thetaE0_2_2, thetaE0_1_3, thetaE0_2_3;
    
    VectorXd thetaS((Nmax - Nmin + 1) / 2);
    VectorXd thetaE((Nmax - Nmin + 1) / 2);
    VectorXd thetaS2(1);
    VectorXd thetaE2(1);
    VectorXd thetaS3(1);
    VectorXd thetaE3(1);
    
    for (Ns_index; Ns_index < Nmax; Ns_index += 2) {
        
        // This two if loop is checking if any set theta_S or theta_E has crossing 0 degree
        
        thetaS0_1 = theta_S(Ns);
        thetaS0_2 = theta_S(Ne);
        thetaE0_1 = theta_E(Ns);
        thetaE0_2 = theta_E(Ne);
        
        
        if (max(theta_S(Ns), theta_S(Ne)) > 3 * M_PI / 2 && min(theta_S(Ns), theta_S(Ne)) < M_PI / 2) {
            thetaS0_1 = min(theta_S(Ns), theta_S(Ne)) + 2 * M_PI;
            thetaS0_2 = max(theta_S(Ns), theta_S(Ne));
            situation = 1;
        }
        
        if (max(theta_E(Ns), theta_E(Ne)) > 3 * M_PI / 2 && min(theta_E(Ns), theta_E(Ne)) < M_PI / 2) {
            thetaE0_1 = min(theta_E(Ns), theta_E(Ne)) + 2 * M_PI;
            thetaE0_2 = max(theta_E(Ns), theta_E(Ne));
            situation = 2;
        }
        
        
        
        thetaS0 = min(thetaS0_1, thetaS0_2);
        thetaE0 = max(thetaE0_1, thetaE0_2);
        
        //Choose the smaller piece
        if (situation == 1) {
            if ((thetaE0 - (thetaS0 - 2 * M_PI)) > (2 * M_PI - (thetaE0 - (thetaS0 - 2 * M_PI)))) {
                
                thetaS(m) = min(thetaE0_1, thetaE0_2);
                thetaE(m) = max(thetaS0_1, thetaS0_2);
                
            }//End of if loop
            else {
                thetaS(m) = thetaS0;
                thetaE(m) = thetaE0 + 2 * M_PI;
            }
            
        }
        else if (situation == 2) {
            if ((thetaE0 - thetaS0) > (2 * M_PI - (thetaE0 - thetaS0))) {
                
                thetaS(m) = min(thetaE0_1, thetaE0_2);
                thetaE(m) = max(thetaS0_1, thetaS0_2) + 2 * M_PI;
                
            }//End of if loop
            else {
                thetaS(m) = thetaS0;
                thetaE(m) = thetaE0;
            }
        }
        else {
            if ((thetaE0 - thetaS0) > (2 * M_PI - (thetaE0 - thetaS0))) {
                
                thetaS(m) = min(thetaE0_1, thetaE0_2);
                thetaE(m) = max(thetaS0_1, thetaS0_2) + 2 * M_PI;
                
            }//End of if loop
            else {
                thetaS(m) = thetaS0;
                thetaE(m) = thetaE0;
            }
        }// End of else loop
        
        m++;
        Ns += 2;
        Ne += 2;
    }//End of for loop
    
    //PC case started.....
    if (PC_term != 0 || APD != 1) {
        
        for (int i = 0; i < nandtheta_array_size2 / 3; i++) {
            
            zoneNum2(i) = nandtheta_array2[0 + i * 3];
            theta_S2(i) = nandtheta_array2[1 + i * 3];
            theta_E2(i) = nandtheta_array2[2 + i * 3];
            
        } // End of for loop
        Nmin02 = (int)zoneNum2.minCoeff();
        Nmax02 = (int)zoneNum2.maxCoeff();
        
        //Define the "real" Nstart & Nend
        if (ns > Nmin02) {
            Nmin2 = ns;
            delta2 = ns - Nmin02;
        }
        else {
            Nmin2 = Nmin02;
            delta2 = 0;
        }
        
        if (ne < Nmax02) {
            Nmax2 = ne;
        }
        else {
            Nmax2 = Nmax02;
        }
        
        // Make sure it start from the same even/odd order as non off-axis case
        
        if ((ns % 2) != 1) {
            if ((Nmin2 % 2) == 1) {
                Nmin2 = Nmin2 + 1;
                delta2 = delta2 + 1;
            }
        }//End of if loop
        else {
            if ((Nmin2 % 2) != 1) {
                Nmin2 = Nmin2 + 1;
                delta2 = delta2 + 1;
            }
        }// End of else loop
        
        // if Nmax is controled by Nmax0, we have to -1 becasue we don't have the
        // information about the starting and ending angle of that zone.
        
        if (((Nmax2 - Nmin2 + 1) % 2) == 1) {
            if (Nmax2 == Nmax02) {
                Nmax2 = Nmax2 - 1;
            }
            else {
                Nmax2 = Nmax2 + 1;
            }
        }//End of if loop
        
        if (Nmin2 >= Nmax2) {
            printf("\n Nmin2 > Nmax2, no ZP allowed\n");
            system("PAUSE");
            return 0;
        }
        
        
        m2 = 0;
        Ns2 = delta2;
        Ne2 = Ns2 + 1;
        Ns_index2 = Nmin2;
        thetaS2.conservativeResize((Nmax2 - Nmin2 + 1) / 2);
        thetaE2.conservativeResize((Nmax2 - Nmin2 + 1) / 2);
        
        
        
        for (Ns_index2; Ns_index2 < Nmax2; Ns_index2 += 2) {
            
            // This two if loop is checking if any set theta_S or theta_E has crossing 0 degree
            
            thetaS0_1_2 = theta_S2(Ns2);
            thetaS0_2_2 = theta_S2(Ne2);
            thetaE0_1_2 = theta_E2(Ns2);
            thetaE0_2_2 = theta_E2(Ne2);
            
            
            if (max(theta_S2(Ns2), theta_S2(Ne2)) > 3 * M_PI / 2 && min(theta_S2(Ns2), theta_S2(Ne2)) < M_PI / 2) {
                thetaS0_1_2 = min(theta_S2(Ns2), theta_S2(Ne2)) + 2 * M_PI;
                thetaS0_2_2 = max(theta_S2(Ns2), theta_S2(Ne2));
                situation2 = 1;
            }
            
            if (max(theta_E2(Ns2), theta_E2(Ne2)) > 3 * M_PI / 2 && min(theta_E2(Ns2), theta_E2(Ne2)) < M_PI / 2) {
                thetaE0_1_2 = min(theta_E2(Ns2), theta_E2(Ne2)) + 2 * M_PI;
                thetaE0_2_2 = max(theta_E2(Ns2), theta_E2(Ne2));
                situation2 = 2;
            }
            
            
            
            thetaS02 = min(thetaS0_1_2, thetaS0_2_2);
            thetaE02 = max(thetaE0_1_2, thetaE0_2_2);
            
            //Choose the smaller piece
            if (situation2 == 1) {
                if ((thetaE02 - (thetaS02 - 2 * M_PI)) > (2 * M_PI - (thetaE02 - (thetaS02 - 2 * M_PI)))) {
                    
                    thetaS2(m2) = min(thetaE0_1_2, thetaE0_2_2);
                    thetaE2(m2) = max(thetaS0_1_2, thetaS0_2_2);
                    
                }//End of if loop
                else {
                    thetaS2(m2) = thetaS02;
                    thetaE2(m2) = thetaE02 + 2 * M_PI;
                }
                
            }
            else if (situation2 == 2) {
                if ((thetaE02 - thetaS02) > (2 * M_PI - (thetaE02 - thetaS02))) {
                    
                    thetaS2(m2) = min(thetaE0_1_2, thetaE0_2_2);
                    thetaE2(m2) = max(thetaS0_1_2, thetaS0_2_2) + 2 * M_PI;
                    
                }//End of if loop
                else {
                    thetaS2(m2) = thetaS02;
                    thetaE2(m2) = thetaE02;
                }
            }
            else {
                if ((thetaE02 - thetaS02) > (2 * M_PI - (thetaE02 - thetaS02))) {
                    
                    thetaS2(m2) = min(thetaE0_1_2, thetaE0_2_2);
                    thetaE2(m2) = max(thetaS0_1_2, thetaS0_2_2) + 2 * M_PI;
                    
                }//End of if loop
                else {
                    thetaS2(m2) = thetaS02;
                    thetaE2(m2) = thetaE02;
                }
            }// End of else loop
            
            m2++;
            Ns2 += 2;
            Ne2 += 2;
        }//End of for loop
        
    }// End of if (PC_term != 0)
    
    //Inner case started.....
    if (InnerSigma != 0) {
        
        for (int i = 0; i < nandtheta_array_size3 / 3; i++) {
            
            zoneNum3(i) = nandtheta_array3[0 + i * 3];
            theta_S3(i) = nandtheta_array3[1 + i * 3];
            theta_E3(i) = nandtheta_array3[2 + i * 3];
            
        } // End of for loop
        Nmin03 = (int)zoneNum3.minCoeff();
        Nmax03 = (int)zoneNum3.maxCoeff();
        
        //Define the "real" Nstart & Nend
        if (ns > Nmin03) {
            Nmin3 = ns;
            delta3 = ns - Nmin03;
        }
        else {
            Nmin3 = Nmin03;
            delta3 = 0;
        }
        
        if (ne < Nmax03) {
            Nmax3 = ne;
        }
        else {
            Nmax3 = Nmax03;
        }
        
        // Make sure it start from the same even/odd order as non off-axis case
        
        if ((ns % 2) != 1) {
            if ((Nmin3 % 2) == 1) {
                Nmin3 = Nmin3 + 1;
                delta3 = delta3 + 1;
            }
        }//End of if loop
        else {
            if ((Nmin3 % 2) != 1) {
                Nmin3 = Nmin3 + 1;
                delta3 = delta3 + 1;
            }
        }// End of else loop
        
        // if Nmax is controled by Nmax0, we have to -1 becasue we don't have the
        // information about the starting and ending angle of that zone.
        
        if (((Nmax3 - Nmin3 + 1) % 2) == 1) {
            if (Nmax3 == Nmax03) {
                Nmax3 = Nmax3 - 1;
            }
            else {
                Nmax3 = Nmax3 + 1;
            }
        }//End of if loop
        
        if (Nmin3 >= Nmax3) {
            printf("\n Nmin3 > Nmax3, no ZP allowed\n");
            system("PAUSE");
            return 0;
        }
        
        
        m3 = 0;
        Ns3 = delta3;
        Ne3 = Ns3 + 1;
        Ns_index3 = Nmin3;
        thetaS3.conservativeResize((Nmax3 - Nmin3 + 1) / 2);
        thetaE3.conservativeResize((Nmax3 - Nmin3 + 1) / 2);
        
        
        
        for (Ns_index3; Ns_index3 < Nmax3; Ns_index3 += 2) {
            
            // This two if loop is checking if any set theta_S or theta_E has crossing 0 degree
            
            thetaS0_1_3 = theta_S3(Ns3);
            thetaS0_2_3 = theta_S3(Ne3);
            thetaE0_1_3 = theta_E3(Ns3);
            thetaE0_2_3 = theta_E3(Ne3);
            
            
            if (max(theta_S3(Ns3), theta_S3(Ne3)) > 3 * M_PI / 2 && min(theta_S3(Ns3), theta_S3(Ne3)) < M_PI / 2) {
                thetaS0_1_3 = min(theta_S3(Ns3), theta_S3(Ne3)) + 2 * M_PI;
                thetaS0_2_3 = max(theta_S3(Ns3), theta_S3(Ne3));
                situation3 = 1;
            }
            
            if (max(theta_E3(Ns3), theta_E3(Ne3)) > 3 * M_PI / 2 && min(theta_E3(Ns3), theta_E3(Ne3)) < M_PI / 2) {
                thetaE0_1_3 = min(theta_E3(Ns3), theta_E3(Ne3)) + 2 * M_PI;
                thetaE0_2_3 = max(theta_E3(Ns3), theta_E3(Ne3));
                situation3 = 2;
            }
            
            
            
            thetaS03 = min(thetaS0_1_3, thetaS0_2_3);
            thetaE03 = max(thetaE0_1_3, thetaE0_2_3);
            
            //Choose the smaller piece
            if (situation3 == 1) {
                if ((thetaE03 - (thetaS03 - 2 * M_PI)) > (2 * M_PI - (thetaE03 - (thetaS03 - 2 * M_PI)))) {
                    
                    thetaS3(m3) = min(thetaE0_1_3, thetaE0_2_3);
                    thetaE3(m3) = max(thetaS0_1_3, thetaS0_2_3);
                    
                }//End of if loop
                else {
                    thetaS3(m3) = thetaS03;
                    thetaE3(m3) = thetaE03 + 2 * M_PI;
                }
                
            }
            else if (situation3 == 2) {
                if ((thetaE03 - thetaS03) > (2 * M_PI - (thetaE03 - thetaS03))) {
                    
                    thetaS3(m3) = min(thetaE0_1_3, thetaE0_2_3);
                    thetaE3(m3) = max(thetaS0_1_3, thetaS0_2_3) + 2 * M_PI;
                    
                }//End of if loop
                else {
                    thetaS3(m3) = thetaS03;
                    thetaE3(m3) = thetaE03;
                }
            }
            else {
                if ((thetaE03 - thetaS03) > (2 * M_PI - (thetaE03 - thetaS03))) {
                    
                    thetaS3(m3) = min(thetaE0_1_3, thetaE0_2_3);
                    thetaE3(m3) = max(thetaS0_1_3, thetaS0_2_3) + 2 * M_PI;
                    
                }//End of if loop
                else {
                    thetaS3(m3) = thetaS03;
                    thetaE3(m3) = thetaE03;
                }
            }// End of else loop
            
            m3++;
            Ns3 += 2;
            Ne3 += 2;
        }//End of for loop
        
    }// End of if (InnerSigma != 0)
    
    //Gather the information for each loop
    //HWcutomZP version 2: varyintg thetaS & thetaE at every zone
    // secant solve parameters:\
    
    int MAXITER = 20;
    double TOLX = 0.001;
    
    int n, test0;
    double rs1, re1, rc1, rs2, re2, rc2, rs3, re3, rc3, dr, x0;
    double delta_APD;
    
    int numZones = (Nmax - Nmin + 1) / 2;
    int totalslice = 0;
    
    //Just for initial value
    VectorXf ReproduceArc(1);
    //VectorXd data(2);
    //int index = 0;
    
    
    
    // Fit the ZP with arc...... based on the algorism in arcfit3 & use secantsovle to
    // get the information of ZP like rs, re, rc ,dr from it.
    
    // Using 2 degree as starting point:\
    
    double totalarcNum = 360 / 2;
    int index = 0;
    int index2 = 0;
    int	PC_on = 0;
    int APD_on = 0;
    
    VectorXi zoneNumber(numZones);
    VectorXd dosage(numZones);
    VectorXd Radius(numZones);
    double zone_Radius;
    double clock;
    int h = 0;
    float R;
    double max_dose;
    double min_dose;
    double IoP2;
    
    IoP2 = IoP;
    
    if (NoP == IoP) {
        //If you want the last part (2nd piece for NoP = 2)
        IoP = 0;
    }
    
    for (int k = 0; k < numZones; k++) {
        
        //These terms need global definition
        //printf("computing zone #%d out of %d zones\n", k+1, numZones);
        // If we are using the WRV file, don't count the curl at GDSII generation function
        if (k % 5 == 0) {
            if (curl_on == 1) {
                sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                        zpID, ((float)(k+1) / ((float)numZones)) / 2);
                system(stringBuffer);
                //printf("Counting p0.5: %0.3f \n", ((float)(k + 1) / ((float)numZones)) / 2);
            }
        }
        
        
        
        /*//Determine the part is generated.
         //For multiple pattering purpose
         if (NoP == 1 && IoP == 1){
         //Continue, nothing happen
         }
         else if (fmod(k+1, NoP) !=  IoP){
         // Has to jump to the next part is interested!
         // Add the difference to it
         if(k == 0){
         //Lock down the first number.
         k = IoP2 - 1;
         }
         else{
         int Rem = fmod(k+1, NoP);
         if (Rem == 0){
         Rem = NoP;
         }
         k = k + (int) (IoP + NoP - Rem);
         }
         }
         */
        
        if (fmod(k + 1, NoP) != IoP) {
            //Bypass this situation
        }
        else {
            
            int arcNum = 0;
            
            VectorXf eachzoner(1);
            VectorXf eachzonedr(1);
            VectorXf eachzonetheta(1);
            VectorXf eachzonedtheta(1);
            VectorXf eachzonexc(1);
            VectorXf eachzoneyc(1);
            
            // Add a if here to determine whether use radial or pin hole seive
            
            int zoneCt = k;
            int no_PC_region;
            
            //Start the iteration:\
            
            double theta = thetaS(k);
            double theta_End = thetaE(k);
            
            double theta_temp;
            double theta_End_temp;
            
            // For SHARP fix axis only
            // Make sure it covers +y axis
            
            if ((theta - M_PI / 2)*(theta_End - M_PI / 2) > 0) {
                //Stick with the orientation of other zoneplates.
                //Max becomes the starting point, min becomes the end.
                
                double theta_temp = theta;
                double theta_End_temp = theta_End;
                
                theta = theta_End_temp - 2 * M_PI;
                theta_End = theta_temp;
                
                
            }
            
            double theta2;
            double theta_End2;
            double theta3;
            double theta_End3;
            double theta_write_portion, theta_write_portion1, theta_write_portion2;
            double theta_empty_portion, theta_empty_portion1, theta_empty_portion2;
            double endp = (1 / totalarcNum) * 2 * M_PI;
            double endp0 = 0;
            double lastpiece = 0;
            double x1, y1,
            x2, y2,
            xm, ym;
            int empty_Num, empty_Num1, empty_Num2;
            double theta_APDS = 0;
            double theta_APDE = 0;
            int APD_index = 0;
            int APD_index2 = 0;
            double theta_window;
            double reverse;
            
            
            // Compute zone number::\
            
            n = Nmin + k * 2;
            
            if (State == 1) {
                // Do nothinig.... original code
            }
            else if (State == 2) {
                // PC offaxis disk
                if (n >= Nmin2 && n < Nmax2) {
                    
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                }
            }
            else if (State == 3) {
                // PC offaxis annular
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                }
                if (n >= Nmin3 && n < Nmax3) {
                    theta3 = thetaS3(index2);
                    theta_End3 = thetaE3(index2);
                    index2++;
                }
            }
            else if (State == 4) {
                // PC APD offaxis disk
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                    theta_write_portion = (theta_End2 - theta2)*APD;
                    theta_empty_portion = (theta_End2 - theta2)*(1 - APD);
                    
                    // PC APD offaxis disk
                    x0 = sqrt(n * lambda * (p*q) / (p + q));
                    rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(n + 1 * lambda * (p*q) / (p + q));
                    re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    rc1 = (rs1 + re1) / 2;
                    
                    //Estimate the min angle for seperation
                    double theta_min = 200 / (rc1 * 1000);
                    empty_Num = (int)floor(theta_empty_portion / theta_min);
                    delta_APD = (theta_empty_portion - empty_Num* theta_min);
                }
            }
            else if (State == 5) {
                // PC APD offaxis annular
                
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                    
                    if (n >= Nmin3 && n < Nmax3) {
                        theta3 = thetaS3(index2);
                        theta_End3 = thetaE3(index2);
                        index2++;
                        // In annular region
                        
                        theta_write_portion1 = (theta3 - theta2)*APD;
                        theta_empty_portion1 = (theta3 - theta2)*(1 - APD);
                        
                        theta_write_portion2 = (theta_End2 - theta_End3)*APD;
                        theta_empty_portion2 = (theta_End2 - theta_End3)*(1 - APD);
                        
                        double theta_min = 200 / (rc1 * 1000);
                        empty_Num1 = (int)floor(theta_empty_portion1 / theta_min);
                        empty_Num2 = (int)floor(theta_empty_portion2 / theta_min);
                    }
                    else {
                        // Haven't hit inner cirlce region
                        theta_write_portion = (theta_End2 - theta2)*APD;
                        theta_empty_portion = (theta_End2 - theta2)*(1 - APD);
                        
                        // PC APD offaxis disk
                        x0 = sqrt(n * lambda * (p*q) / (p + q));
                        rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                          PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(n + 1 * lambda * (p*q) / (p + q));
                        re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                          PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        rc1 = (rs1 + re1) / 2;
                        
                        //Estimate the min angle for seperation
                        double theta_min = 200 / (rc1 * 1000);
                        empty_Num = (int)floor(theta_empty_portion / theta_min);
                    }
                }
            }
            else if (State == 6) {
                // APD offaxis disk
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                    theta_write_portion = (theta_End2 - theta2)*APD;
                    theta_empty_portion = (theta_End2 - theta2)*(1 - APD);
                    
                    // APD offaxis disk
                    x0 = sqrt(n * lambda * (p*q) / (p + q));
                    rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER,
                                      TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(n + 1 * lambda * (p*q) / (p + q));
                    re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    rc1 = (rs1 + re1) / 2;
                    
                    //Estimate the min angle for seperation
                    double theta_min = 200 / (rc1 * 1000);
                    empty_Num = (int)floor(theta_empty_portion / theta_min);
                    delta_APD = (theta_empty_portion - empty_Num* theta_min);
                }
            }
            
            else if (State == 7) {
                // APD offaxis annular
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                    
                    if (n >= Nmin3 && n < Nmax3) {
                        theta3 = thetaS3(index2);
                        theta_End3 = thetaE3(index2);
                        index2++;
                        // In annular region
                        
                        theta_write_portion1 = (theta3 - theta2)*APD;
                        theta_empty_portion1 = (theta3 - theta2)*(1 - APD);
                        
                        theta_write_portion2 = (theta_End2 - theta_End3)*APD;
                        theta_empty_portion2 = (theta_End2 - theta_End3)*(1 - APD);
                        
                        double theta_min = 200 / (rc1 * 1000);
                        empty_Num1 = (int)floor(theta_empty_portion1 / theta_min);
                        empty_Num2 = (int)floor(theta_empty_portion2 / theta_min);
                    }
                    else {
                        // Haven't hit inner cirlce region
                        theta_write_portion = (theta_End2 - theta2)*APD;
                        theta_empty_portion = (theta_End2 - theta2)*(1 - APD);
                        
                        // PC APD offaxis disk
                        x0 = sqrt(n * lambda * (p*q) / (p + q));
                        rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                          PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(n + 1 * lambda * (p*q) / (p + q));
                        re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                          PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        rc1 = (rs1 + re1) / 2;
                        
                        //Estimate the min angle for seperation
                        double theta_min = 200 / (rc1 * 1000);
                        empty_Num = (int)floor(theta_empty_portion / theta_min);
                    }
                }
            }
            else if (State == 8) {
                //Dark field disk
                // Only need to know the starting and ending point
                if (n >= Nmin2 && n < Nmax2) {
                    
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                }
            }
            else if (State == 9) {
                // Dark field annular
                if (n >= Nmin2 && n < Nmax2) {
                    theta2 = thetaS2(index);
                    theta_End2 = thetaE2(index);
                    index++;
                }
                if (n >= Nmin3 && n < Nmax3) {
                    theta3 = thetaS3(index2);
                    theta_End3 = thetaE3(index2);
                    index2++;
                }
            }
            
            //For state = 10 only APD on offaxis ZP
            double jump2 = 1;
            
            while (theta < theta_End - 0.00001) {
                
                PC_on = 0;
                APD_on = 0;
                
                if (State == 1) {
                    // offaxis do nothing......
                }
                else if (State == 2) {
                    // PC offaxis disk
                    if (n >= Nmin2 && n < Nmax2) {
                        if (theta >= theta2 && theta < theta_End2) {
                            PC_on = 1;
                            no_PC_region = 0;
                            lastpiece = 0;
                        }
                    }
                }
                else if (State == 3) {
                    // PC offaxis annular
                    if (n >= Nmin2 && n < Nmax2) {
                        // If it is annular region
                        if (n >= Nmin3 && n < Nmax3) {
                            if (theta >= theta2 && theta < theta3) {
                                PC_on = 1;
                                no_PC_region = 0;
                                lastpiece = 0;
                            }
                            if (theta >= theta_End3 && theta < theta_End2) {
                                PC_on = 1;
                                no_PC_region = 0;
                                lastpiece = 0;
                            }
                        }
                        else {// Hasn't hit the inner circle
                            if (theta >= theta2 && theta < theta_End2) {
                                PC_on = 1;
                                no_PC_region = 0;
                                lastpiece = 0;
                            }
                        }
                    }
                }
                else if (State == 4) {
                    // PC APD offaxis disk
                    if (n >= Nmin2 && n < Nmax2) {
                        if (theta >= theta2 && theta < theta_End2) {
                            if (APD_index == 0) {
                                theta_APDS = theta2;
                                theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                            }
                            else {
                                theta_APDS = theta2 + APD_index* ((theta_empty_portion + theta_write_portion) / empty_Num);
                                theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                            }
                            if (theta >= theta_APDS && theta< theta_APDE) {
                                PC_on = 1;
                                no_PC_region = 0;
                                lastpiece = 0;
                            }
                            else {
                                if (theta >= theta2 && theta < theta_End2) {
                                    // Check if it is the last piece
                                    if (APD_index == empty_Num) {
                                        theta = theta_End2;
                                    }
                                    else {
                                        // hitting the empty region jump to next one
                                        if (theta < theta_APDS) {
                                            theta = theta_APDS;
                                            PC_on = 1;
                                            no_PC_region = 0;
                                            lastpiece = 0;
                                        }
                                        // hitting the empty region jump to next one
                                        if (theta >= theta_APDE) {
                                            theta = theta_APDE + theta_empty_portion / empty_Num;
                                            PC_on = 1;
                                            no_PC_region = 0;
                                        }
                                    }
                                }
                            }
                            // if it is the last piece, hitting the edge of theta_End2
                            // turn off PC_term
                            
                            if (theta >= theta_End2) {
                                PC_on = 0;
                                no_PC_region = 1;
                            }
                            
                            APD_index++;
                        }
                    }
                }
                else if (State == 5) {
                    //PC APD offaxis annular
                    if (n >= Nmin2 && n < Nmax2) {
                        if (n >= Nmin3 && n < Nmax3) {
                            // In the 1st annular region
                            if (theta >= theta2 && theta < theta3) {
                                if (APD_index == 0) {
                                    theta_APDS = theta2;
                                    theta_APDE = theta_APDS + theta_write_portion1 / empty_Num1;
                                }
                                else {
                                    theta_APDS = theta2 + APD_index* ((theta_empty_portion1 + theta_write_portion1) / empty_Num1);
                                    theta_APDE = theta_APDS + theta_write_portion1 / empty_Num1;
                                }
                                lastpiece = 0;
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 1;
                                    no_PC_region = 0;
                                }
                                else {
                                    // hitting the empty region jump to next one
                                    if (theta < theta_APDS) {
                                        theta = theta_APDS;
                                        PC_on = 1;
                                        no_PC_region = 0;
                                    }
                                    // hitting the empty region jump to next one
                                    if (theta >= theta_APDE) {
                                        theta = theta_APDE + theta_empty_portion1 / empty_Num1;
                                        PC_on = 1;
                                        no_PC_region = 0;
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta3) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index++;
                            }
                            
                            // In the 2nd annular region
                            if (theta >= theta_End3 && theta < theta_End2) {
                                if (APD_index2 == 0) {
                                    theta_APDS = theta_End3;
                                    theta_APDE = theta_APDS + theta_write_portion2 / empty_Num2;
                                }
                                else {
                                    theta_APDS = theta_End3 + APD_index2* ((theta_empty_portion2 + theta_write_portion2) / empty_Num2);
                                    theta_APDE = theta_APDS + theta_write_portion2 / empty_Num2;
                                }
                                
                                lastpiece = 0;
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 1;
                                    no_PC_region = 0;
                                    lastpiece = 0;
                                }
                                else {
                                    // hitting the empty region jump to next one
                                    if (theta < theta_APDS) {
                                        theta = theta_APDS;
                                        PC_on = 1;
                                        no_PC_region = 0;
                                        lastpiece = 0;
                                    }
                                    // hitting the empty region jump to next one
                                    if (theta >= theta_APDE) {
                                        theta = theta_APDE + theta_empty_portion2 / empty_Num2;
                                        PC_on = 1;
                                        no_PC_region = 0;
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta_End2) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index2++;
                            }
                        }
                        else {
                            if (theta >= theta2 && theta < theta_End2) {
                                if (APD_index == 0) {
                                    theta_APDS = theta2;
                                    theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                                }
                                else {
                                    theta_APDS = theta2 + APD_index* ((theta_empty_portion + theta_write_portion) / empty_Num);
                                    theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                                }
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 1;
                                    no_PC_region = 0;
                                    lastpiece = 0;
                                }
                                else {
                                    if (theta >= theta2 && theta < theta_End2) {
                                        // Check if it is the last piece
                                        if (APD_index == empty_Num) {
                                            theta = theta_End2;
                                        }
                                        else {
                                            // hitting the empty region jump to next one
                                            if (theta < theta_APDS) {
                                                theta = theta_APDS;
                                                PC_on = 1;
                                                no_PC_region = 0;
                                                lastpiece = 0;
                                            }
                                            // hitting the empty region jump to next one
                                            if (theta >= theta_APDE) {
                                                theta = theta_APDE + theta_empty_portion / empty_Num;
                                                PC_on = 1;
                                                no_PC_region = 0;
                                            }
                                        }
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta_End2) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index++;
                            }
                        }// End of else loop
                    }// if (n >= Nmin2 && n < Nmax2)
                    
                }// End of state 5
                else if (State == 6) {
                    // APD offaxis disk
                    // PC_on == 0
                    if (n >= Nmin2 && n < Nmax2) {
                        if (theta >= theta2 && theta < theta_End2) {
                            if (APD_index == 0) {
                                theta_APDS = theta2;
                                theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                            }
                            else {
                                theta_APDS = theta2 + APD_index* ((theta_empty_portion + theta_write_portion) / empty_Num);
                                theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                            }
                            if (theta >= theta_APDS && theta< theta_APDE) {
                                PC_on = 0;
                                no_PC_region = 1;
                                lastpiece = 0;
                            }
                            else {
                                if (theta >= theta2 && theta < theta_End2) {
                                    // Check if it is the last piece
                                    if (APD_index == empty_Num) {
                                        theta = theta_End2;
                                    }
                                    else {
                                        // hitting the empty region jump to next one
                                        if (theta < theta_APDS) {
                                            theta = theta_APDS;
                                            PC_on = 0;
                                            no_PC_region = 1;
                                            lastpiece = 0;
                                        }
                                        // hitting the empty region jump to next one
                                        if (theta >= theta_APDE) {
                                            theta = theta_APDE + theta_empty_portion / empty_Num;
                                            PC_on = 0;
                                            no_PC_region = 1;
                                        }
                                    }
                                }
                            }
                            // if it is the last piece, hitting the edge of theta_End2
                            // turn off PC_term
                            
                            if (theta >= theta_End2) {
                                PC_on = 0;
                                no_PC_region = 1;
                            }
                            
                            APD_index++;
                        }
                    }
                }// End of state 6
                else if (State == 7) {
                    // APD offaxis annular
                    // PC_on == 0
                    if (n >= Nmin2 && n < Nmax2) {
                        if (n >= Nmin3 && n < Nmax3) {
                            // In the 1st annular region
                            if (theta >= theta2 && theta < theta3) {
                                if (APD_index == 0) {
                                    theta_APDS = theta2;
                                    theta_APDE = theta_APDS + theta_write_portion1 / empty_Num1;
                                }
                                else {
                                    theta_APDS = theta2 + APD_index* ((theta_empty_portion1 + theta_write_portion1) / empty_Num1);
                                    theta_APDE = theta_APDS + theta_write_portion1 / empty_Num1;
                                }
                                lastpiece = 0;
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 1;
                                    no_PC_region = 0;
                                }
                                else {
                                    // hitting the empty region jump to next one
                                    if (theta < theta_APDS) {
                                        theta = theta_APDS;
                                        PC_on = 0;
                                        no_PC_region = 1;
                                    }
                                    // hitting the empty region jump to next one
                                    if (theta >= theta_APDE) {
                                        theta = theta_APDE + theta_empty_portion1 / empty_Num1;
                                        PC_on = 0;
                                        no_PC_region = 1;
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta3) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index++;
                            }
                            
                            // In the 2nd annular region
                            if (theta >= theta_End3 && theta < theta_End2) {
                                if (APD_index2 == 0) {
                                    theta_APDS = theta_End3;
                                    theta_APDE = theta_APDS + theta_write_portion2 / empty_Num2;
                                }
                                else {
                                    theta_APDS = theta_End3 + APD_index2* ((theta_empty_portion2 + theta_write_portion2) / empty_Num2);
                                    theta_APDE = theta_APDS + theta_write_portion2 / empty_Num2;
                                }
                                
                                lastpiece = 0;
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                    lastpiece = 0;
                                }
                                else {
                                    // hitting the empty region jump to next one
                                    if (theta < theta_APDS) {
                                        theta = theta_APDS;
                                        PC_on = 0;
                                        no_PC_region = 1;
                                        lastpiece = 0;
                                    }
                                    // hitting the empty region jump to next one
                                    if (theta >= theta_APDE) {
                                        theta = theta_APDE + theta_empty_portion2 / empty_Num2;
                                        PC_on = 0;
                                        no_PC_region = 1;
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta_End2) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index2++;
                            }
                        }
                        else {
                            if (theta >= theta2 && theta < theta_End2) {
                                if (APD_index == 0) {
                                    theta_APDS = theta2;
                                    theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                                }
                                else {
                                    theta_APDS = theta2 + APD_index* ((theta_empty_portion + theta_write_portion) / empty_Num);
                                    theta_APDE = theta_APDS + theta_write_portion / empty_Num;
                                }
                                
                                if (theta >= theta_APDS && theta< theta_APDE) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                    lastpiece = 0;
                                }
                                else {
                                    if (theta >= theta2 && theta < theta_End2) {
                                        // Check if it is the last piece
                                        if (APD_index == empty_Num) {
                                            theta = theta_End2;
                                        }
                                        else {
                                            // hitting the empty region jump to next one
                                            if (theta < theta_APDS) {
                                                theta = theta_APDS;
                                                PC_on = 0;
                                                no_PC_region = 1;
                                                lastpiece = 0;
                                            }
                                            // hitting the empty region jump to next one
                                            if (theta >= theta_APDE) {
                                                theta = theta_APDE + theta_empty_portion / empty_Num;
                                                PC_on = 0;
                                                no_PC_region = 1;
                                            }
                                        }
                                    }
                                }
                                // if it is the last piece, hitting the edge of theta_End2
                                // turn of PC_term
                                
                                if (theta >= theta_End2) {
                                    PC_on = 0;
                                    no_PC_region = 1;
                                }
                                
                                APD_index++;
                            }
                        }// End of else loop
                    }// if (n >= Nmin2 && n < Nmax2)
                    
                }// End of state 7
                
                else if (State == 8) {
                    // Darkfield offaxis disk
                    if (n >= Nmin2 && n < Nmax2) {
                        if (theta >= theta2 && theta < theta_End2) {
                            PC_on = 0;
                            no_PC_region = 1;
                            lastpiece = 0;
                            theta = theta_End2;
                        }
                    }
                    
                }// End of state 8
                else if (State == 9) {
                    // Darkfield offaxis annular
                    if (n >= Nmin2 && n < Nmax2) {
                        // If it is annular region
                        if (n >= Nmin3 && n < Nmax3) {
                            if (theta >= theta2 && theta < theta3) {
                                PC_on = 0;
                                no_PC_region = 1;
                                lastpiece = 0;
                                theta = theta3;
                            }
                            if (theta >= theta_End3 && theta < theta_End2) {
                                PC_on = 0;
                                no_PC_region = 1;
                                lastpiece = 0;
                                theta = theta_End2;
                            }
                        }
                        else {// Hasn't hit the inner circle
                            if (theta >= theta2 && theta < theta_End2) {
                                PC_on = 0;
                                no_PC_region = 1;
                                lastpiece = 0;
                                theta = theta_End2;
                            }
                        }
                    }
                    
                }// End of state 9
                else if (State == 10) {
                    //APD window ZP,  define theta and theta_end, and theta_next for each segment:
                    if (jump2 == 1) {
                        x0 = sqrt(n * lambda * (p*q) / (p + q));
                        rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER,
                                          TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(n + 1 * lambda * (p*q) / (p + q));
                        re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                          PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        rc1 = (rs1 + re1) / 2;
                        
                        //Estimate the min angle for seperation, now the min piece is dr by dr
                        double theta_dr = (re1 - rs1) / (rc1);
                        
                        theta_window = theta_dr;
                        
                        double FirstTime = 0;
                        
                        reverse = 0;
                        double theta_window2 = 0;
                        double APD_Wvalue;
                        double miss = 0;
                        
                        while (1) {
                            
                            //Calculate its apodization value based on the cooridnates
                            //APD_Wvalue = transmission
                            //always print the bigger one to avoid the problem
                            if (FirstTime == 0) {
                                APD_Wvalue = APD_WindowF(rc1, theta + theta_window, orders, numA, n, p, lambda, q, R0, PC_term, PC_on, offaxis_aberration, CRA, NA_C, APD_window);
                                
                                //Debug
                                // printf("%f DC\n", APD_Wvalue);
                                
                                
                                //Check the duty cycle
                                if (APD_Wvalue >= 0.5) {
                                    // if it is above 0.5, add one gap, the maximum gap will be dr which will give you 50%
                                    reverse = 1;
                                    theta_window2 = theta_window + theta_dr*APD_Wvalue / (1 - APD_Wvalue);
                                }
                                else {
                                    // Smaller than 50% duty cycle: add one block, to make sure that the block will be larger than dr, avoiding the fabrication limit
                                    reverse = 0;
                                    theta_window2 = theta_window + theta_dr*(1 - APD_Wvalue) / APD_Wvalue;
                                }
                                
                                
                                FirstTime = 1;
                            }
                            
                            //Obtain the new APD_window value and make a comparison
                            double APD_Wvalue2 = APD_WindowF(rc1, theta + theta_window2, orders, numA, n, p, lambda, q, R0, PC_term, PC_on, offaxis_aberration, CRA, NA_C, APD_window);
                            double  APD_WvalueM = (APD_Wvalue2 + APD_Wvalue) / 2;
                            
                            //If satisfy the condition:
                            if (abs(APD_Wvalue2 - APD_WvalueM) + abs(APD_Wvalue - APD_WvalueM) <= 0.01) {
                                if (reverse == 0) {
                                    theta_End2 = theta + theta_window + theta_dr*(1 - APD_WvalueM) / APD_WvalueM;
                                    break;
                                }
                                else {
                                    //if reverse == 1
                                    if (miss == 0) {
                                        theta_End2 = theta + theta_window + theta_dr*APD_WvalueM / (1 - APD_WvalueM);
                                        theta_window = theta_dr*APD_WvalueM / (1 - APD_WvalueM);
                                        break;
                                    }
                                    else {
                                        theta_End2 = theta + theta_window2;
                                        theta_window = theta_window2 - theta_window;
                                        break;
                                    }
                                    //theta = theta + theta_window;
                                }
                            }
                            else {
                                //If it doesn't match the condition, make it closer and check it again.
                                theta_window2 = 0.99* theta_window2;
                                //If you miss, then the new window should be theta_window2, instead of the calculation result;
                                miss = 1;
                            }
                            
                        } // End of while (1)
                        jump2 = 0;
                        lastpiece = 0;
                    }//End of if (jump2 == 1)
                }// End of state 10
                
                
                
                //Debug
                //double Debug = theta_End2 - theta;
                //printf("%f radians\n", Debug);
                
                // Define parameters(test):\
                
                double meettol = 0;
                int test = 0;
                double Stop = 0;
                
                VectorXd Rarc(2);
                VectorXd xc(2);
                VectorXd yc(2);
                VectorXd thetada(2);
                VectorXd thetada0(2);
                VectorXd thetaa0(2);
                
                // Define the starting points
                
                double ns0 = n;
                double ne0 = n + 1;
                double jump = 0;
                
                //Compute rs1 & re1:\
                // If state == 1 && with Zernike aberration
                // You have to redefine the cooridnate system: x0 and theta
                
                x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                rs1 = secantSolve(x0, theta, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                  PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                re1 = secantSolve(x0, theta, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                  PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                
                rc1 = (rs1 + re1) / 2;
                
                x1 = 1000 * rc1* cos(theta);
                y1 = 1000 * rc1* sin(theta);
                
                while (1) {
                    
                    if (State == 1) {
                        // offaxis do nothing......
                    }
                    else if (State == 2) {
                        // PC offaxis disk
                        if (n >= Nmin2 && n < Nmax2) {
                            // theta+endp get into PC region
                            if (theta < theta2 && (theta + endp) >= theta2) {
                                endp = theta2 - theta;
                                Stop = 1;
                                lastpiece = 0;
                            }
                            
                            if (theta >= theta2 && theta < theta_End2) {
                                if (theta + endp >= theta_End2) {
                                    endp = theta_End2 - theta;
                                    Stop = 1;
                                    lastpiece = 0;
                                }
                            }
                        }
                    }
                    else if (State == 3) {
                        // PC offaxis annular
                        if (n >= Nmin2 && n < Nmax2) {
                            if (n >= Nmin3 && n < Nmax3) {
                                if (theta < theta2) {
                                    if ((theta + endp) >= theta2 && (theta + endp)< theta3) {
                                        endp = theta2 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                if (theta >= theta2 && theta < theta3) {
                                    if ((theta + endp) >= theta3) {
                                        endp = theta3 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                if (theta >= theta3 && theta < theta_End3) {
                                    if ((theta + endp) >= theta_End3) {
                                        endp = theta_End3 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                if (theta >= theta_End3 && theta < theta_End2) {
                                    if ((theta + endp) >= theta_End2) {
                                        endp = theta_End2 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                            }
                            else {
                                //Haven't hit the inner circle yet
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    endp = theta2 - theta;
                                    Stop = 1;
                                    lastpiece = 0;
                                }
                                
                                if (theta >= theta2 && theta < theta_End2) {
                                    if (theta + endp >= theta_End2) {
                                        endp = theta_End2 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                            }
                        }
                    }
                    else if (State == 4) {
                        // PC APD offaxis disk
                        if (n >= Nmin2 && n < Nmax2) {
                            if (theta < theta2 && (theta + endp) >= theta2) {
                                Stop = 1;
                                endp = theta2 - theta;
                                lastpiece = 0;
                            }
                            if (theta >= theta2 && theta < theta_End2) {
                                if (theta >= theta_APDS && theta <theta_APDE) {
                                    if (theta + endp >= theta_APDE) {
                                        endp = theta_APDE - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                            }
                            if (theta > theta_End2 && (theta + endp) >= theta_End) {
                                Stop = 1;
                                endp = theta_End - theta;
                                lastpiece = 0;
                            }
                        }
                    }
                    else if (State == 5) {
                        if (n >= Nmin2 && n < Nmax2) {
                            if (n >= Nmin3 && n < Nmax3) {
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    endp = theta2 - theta;
                                    Stop = 1;
                                    lastpiece = 0;
                                }
                                // In the 1st region
                                if (theta >= theta2 && theta< theta3) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if ((theta + endp) >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                                
                                if (theta >= theta3 && theta < theta_End3) {
                                    if ((theta + endp) >= theta_End3) {
                                        endp = theta_End3 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                if (theta >= theta_End3 && theta < theta_End2) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if ((theta + endp) >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                            }
                            else {
                                // Not hit the inner region yet
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    Stop = 1;
                                    endp = theta2 - theta;
                                    lastpiece = 0;
                                }
                                if (theta >= theta2 && theta < theta_End2) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if (theta + endp >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                                if (theta > theta_End2 && (theta + endp) >= theta_End) {
                                    Stop = 1;
                                    endp = theta_End - theta;
                                    lastpiece = 0;
                                }
                            }
                        }
                    }
                    else if (State == 6) {
                        // APD offaxis disk
                        if (n >= Nmin2 && n < Nmax2) {
                            if (theta < theta2 && (theta + endp) >= theta2) {
                                Stop = 1;
                                endp = theta2 - theta;
                                lastpiece = 0;
                            }
                            if (theta >= theta2 && theta < theta_End2) {
                                if (theta >= theta_APDS && theta <theta_APDE) {
                                    if (theta + endp >= theta_APDE) {
                                        endp = theta_APDE - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                            }
                            if (theta > theta_End2 && (theta + endp) >= theta_End) {
                                Stop = 1;
                                endp = theta_End - theta;
                                lastpiece = 0;
                            }
                        }
                    }
                    else if (State == 7) {
                        if (n >= Nmin2 && n < Nmax2) {
                            if (n >= Nmin3 && n < Nmax3) {
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    endp = theta2 - theta;
                                    Stop = 1;
                                    lastpiece = 0;
                                }
                                // In the 1st region
                                if (theta >= theta2 && theta< theta3) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if ((theta + endp) >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                                
                                if (theta >= theta3 && theta < theta_End3) {
                                    if ((theta + endp) >= theta_End3) {
                                        endp = theta_End3 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                if (theta >= theta_End3 && theta < theta_End2) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if ((theta + endp) >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                            }
                            else {
                                // Not hit the inner region yet
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    Stop = 1;
                                    endp = theta2 - theta;
                                    lastpiece = 0;
                                }
                                if (theta >= theta2 && theta < theta_End2) {
                                    if (theta >= theta_APDS && theta <theta_APDE) {
                                        if (theta + endp >= theta_APDE) {
                                            endp = theta_APDE - theta;
                                            Stop = 1;
                                            lastpiece = 0;
                                        }
                                    }
                                }
                                if (theta > theta_End2 && (theta + endp) >= theta_End) {
                                    Stop = 1;
                                    endp = theta_End - theta;
                                    lastpiece = 0;
                                }
                            }
                        }
                    }
                    else if (State == 8) {
                        // Darkfield offaxis disk
                        if (n >= Nmin2 && n < Nmax2) {
                            // theta+endp get into Darkfield region
                            if (theta < theta2 && (theta + endp) >= theta2) {
                                endp = theta2 - theta;
                                Stop = 1;
                                lastpiece = 0;
                            }
                        }
                    }
                    else if (State == 9) {
                        //Darkfield annular
                        if (n >= Nmin2 && n < Nmax2) {
                            if (n >= Nmin3 && n < Nmax3) {
                                if (theta < theta2) {
                                    if ((theta + endp) >= theta2 && (theta + endp)< theta3) {
                                        endp = theta2 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                                
                                if (theta >= theta3 && theta < theta_End3) {
                                    if ((theta + endp) >= theta_End3) {
                                        endp = theta_End3 - theta;
                                        Stop = 1;
                                        lastpiece = 0;
                                    }
                                }
                            }
                            else {
                                //Haven't hit the inner circle yet
                                if (theta < theta2 && (theta + endp) >= theta2) {
                                    endp = theta2 - theta;
                                    Stop = 1;
                                    lastpiece = 0;
                                }
                                
                            }
                        }
                    }
                    else if (State == 10) {
                        // In case it get to the gap region
                        if (reverse == 0) {
                            if ((theta + endp) - (theta + theta_window) >= -0.00001) {
                                endp = theta_window;
                                lastpiece = 1;
                                jump = 1;
                            }
                        }
                        else {
                            //if ((theta + endp) - theta_End2 >= -0.00001){
                            if ((theta + endp) - (theta + theta_window) >= -0.00001) {
                                //Has to get rid of the case that only small fracture is not ready,
                                //You don't want to cut a super small piece to fit it!
                                //endp = theta_End2 - theta;
                                endp = theta_window;
                                lastpiece = 1;
                                jump = 1;
                                //lastpiece = 1;
                            }
                        }
                    }
                    // Define the ending point:\
                    // If state == 1 && with Zernike aberration
                    // You have to redefine the coordinate system: x0, theta+endp
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs2 = secantSolve(x0, theta + endp, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re2 = secantSolve(x0, theta + endp, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc2 = (rs2 + re2) / 2;
                    
                    x2 = 1000 * rc2* cos(theta + endp);
                    y2 = 1000 * rc2* sin(theta + endp);
                    
                    // Define the centre point:\
                    // If state == 1 && with Zernike aberration
                    // You have to redefine the coordinate system: x0, theta+endp
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs3 = secantSolve(x0, theta + (endp / 2), orders, numA, ns0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re3 = secantSolve(x0, theta + (endp / 2), orders, numA, ne0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                      PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc3 = (rs3 + re3) / 2;
                    
                    xm = 1000 * rc3* cos(theta + (endp / 2));
                    ym = 1000 * rc3* sin(theta + (endp / 2));
                    
                    //Find the arc with its r, xc, yc, theta, and dtheta based on these three points:\
                    
                    double xsm = (x1 + xm) / 2;
                    double ysm = (y1 + ym) / 2;
                    
                    double xme = (xm + x2) / 2;
                    double yme = (ym + y2) / 2;
                    
                    double slopeSM = (ym - y1) / (xm - x1);
                    double slopeME = (y2 - ym) / (x2 - xm);
                    
                    xc(test) = ((slopeSM*xme - slopeME*xsm) + (yme - ysm)*slopeSM*slopeME) / (slopeSM - slopeME);
                    yc(test) = (xsm - xme + ysm*slopeSM - yme*slopeME) / (slopeSM - slopeME);
                    
                    Rarc(test) = sqrt((x1 - xc(test))*(x1 - xc(test)) + (y1 - yc(test))*(y1 - yc(test)));
                    
                    // Define the coordinate based on xc yc:\
                    
                    double x_arc = x1 - xc(test);
                    double y_arc = y1 - yc(test);
                    
                    //theta: Define as the angle based on the xc & yc to the middle point of the arc:\
                    
                    if (x_arc >= 0 && y_arc >= 0) {
                        thetaa0(test) = asin(y_arc / Rarc(test));
                    }
                    else if (x_arc<0 && y_arc >= 0) {
                        thetaa0(test) = M_PI - asin(y_arc / Rarc(test));
                    }
                    else if (x_arc <= 0 && y_arc<0) {
                        thetaa0(test) = M_PI + asin(abs(y_arc / Rarc(test)));
                    }
                    else {
                        thetaa0(test) = 2 * M_PI - asin(abs(y_arc / Rarc(test)));
                    }
                    
                    //dtheta:\
                    
                    double x_arca = x2 - xc(test);
                    double y_arca = y2 - yc(test);
                    
                    if (x_arca >= 0 && y_arca >= 0) {
                        thetada0(test) = asin(y_arca / Rarc(test));
                    }
                    else if (x_arca<0 && y_arca >= 0) {
                        thetada0(test) = M_PI - asin(y_arca / Rarc(test));
                    }
                    else if (x_arca <= 0 && y_arca<0) {
                        thetada0(test) = M_PI + asin(abs(y_arca / Rarc(test)));
                    }
                    else {
                        thetada0(test) = 2 * M_PI - asin(abs(y_arca / Rarc(test)));
                    }
                    //make sure thetaa0 is always bigger than theta0:\
                    
                    
                    if (thetada0(test) > thetaa0(test)) {
                        thetada(test) = thetada0(test) - thetaa0(test);
                    }
                    else {
                        thetada(test) = thetada0(test) + (2 * M_PI - thetaa0(test));
                    }
                    
                    // Reproduce all the points in the arc and by Reproduce, ready for comparison and iteration!:\
                    
                    MatrixXd coord1(2, 11);
                    VectorXd xa(11);
                    VectorXd ya(11);
                    VectorXd th(11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        th(k) = thetada(test)*k / 10;
                        
                        xa(k) = Rarc(test)*cos(thetaa0(test) + th(k)) + xc(test);
                        ya(k) = Rarc(test)*sin(thetaa0(test) + th(k)) + yc(test);
                        
                        coord1(0, k) = xa(k);
                        coord1(1, k) = ya(k);
                    }
                    
                    MatrixXd coord2(2, 11);
                    VectorXd xr(11);
                    VectorXd yr(11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        double th = endp*k / 10;
                        
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        double rsr = secantSolve(x0, theta + th, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                                 PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        double rer = secantSolve(x0, theta + th, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                                 PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        double rcr = (rsr + rer) / 2;
                        
                        xr(k) = 1000 * rcr* cos(theta + th);
                        yr(k) = 1000 * rcr* sin(theta + th);
                        
                        coord2(0, k) = xr(k);
                        coord2(1, k) = yr(k);
                    }
                    
                    // Based on coord1 & coord2, comparing the difference with tolerance to start the iteration:\
                    
                    VectorXd output(11);
                    
                    for (int k = 0; k <11; k++) {
                        output(k) = sqrt((coord1(0, k) - coord2(0, k))*(coord1(0, k) - coord2(0, k)) + (coord1(1, k) - coord2(1, k))*(coord1(1, k) - coord2(1, k)));
                    }
                    
                    double result = output.maxCoeff();
                    double	dth = endp / 10;
                    double rs, re;
                    double drM = 0;
                    
                    
                    for (int j = 0; j < 11; j++) {
                        double thetaj = theta + j*dth;
                        
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        rs = secantSolve(x0, thetaj, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                         PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        re = secantSolve(x0, thetaj, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, R0, 1, MAXITER, TOLX,
                                         PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        drM += (re - rs);
                    }
                    
                    double result2 = 1000 * drM / 11;
                    double condition = result / result2;
                    
                    if (condition > zTol) {
                        if (meettol == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, R0, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            theta = theta + endp0;
                            test0 = (test == 0);
                            if (theta + endp0 > theta_End) {
                                endp = theta_End - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = endp0;
                                lastpiece = 0;
                            }
                            
                            break;
                        }
                        
                        else {
                            endp = 0.99*endp;
                            lastpiece = 0;
                            test = (test == 0);
                        }
                    }
                    else {
                        if (lastpiece == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, R0, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            test0 = test;
                            
                            if (jump == 1) {
                                //For state = 10, offaxis ZP APD situation only
                                //Jump2 is used to indicate the "while  theta < theta_End" loop to work on the next one
                                theta = theta_End2;
                                jump2 = 1;
                            }
                            else {
                                theta = theta + endp;
                                //theta = theta_End2;
                                //jump2 = 1;
                            }
                            
                            break;
                        }
                        else if (Stop == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, R0, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            test0 = test;
                            theta = theta + endp;
                            break;
                        }
                        else {
                            endp0 = endp;
                            meettol = 1;
                            test = (test == 0);
                            if (theta + 1.01*endp > theta_End) {
                                endp = theta_End - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = 1.01*endp;
                            }
                        }
                    }
                }//end of while 1:\
                
                // Switch back to um: 03042014 by Henry Wang
                
                R = (Rarc(test0) - dr / 2) / 1000;
                eachzoner(arcNum) = (Rarc(test0) - dr / 2) / 1000;
                eachzonedr(arcNum) = dr / 1000;
                eachzonetheta(arcNum) = thetaa0(test0) + (thetada(test0) / 2);
                eachzonedtheta(arcNum) = thetada(test0);
                eachzonexc(arcNum) = xc(test0) / 1000;
                eachzoneyc(arcNum) = yc(test0) / 1000;
                
                arcNum++;
                
                eachzoner.conservativeResize(1 + arcNum);
                eachzonedr.conservativeResize(1 + arcNum);
                eachzonetheta.conservativeResize(1 + arcNum);
                eachzonedtheta.conservativeResize(1 + arcNum);
                eachzonexc.conservativeResize(1 + arcNum);
                eachzoneyc.conservativeResize(1 + arcNum);
                
            }//End of while theta < 2pi;:\
            
            printf("Finished zone with %d arcs.\n", arcNum);
            
            if (file_format == 4 || file_format == 2) {
                
                zoneNumber(h) = n;
                //dosage(h) = (dr/1000+ bias/1000)/(dr/1000+ bias/1000- bias/1000);
                dosage(h) = (dr / 1000) / (dr / 1000 - bias / 1000);
                Radius(h) = R;
                h = h + 1;
                
            }
            
            // Save all the information for each zones
            // We need to save 6 different parameters (rows), and all the information
            // from each zone is saved in columns:\
            
            int totalslice0 = totalslice;
            ReproduceArc.conservativeResize((totalslice + arcNum) * 6);
            
            for (int l = 0; l < arcNum; l++) {
                
                // Need to save the data in the form of "{R, theta, dR, dtheta, Xc, Yc }":\
                // Due to fabrication consideration, we need to add bias term to reduce the zone width from
                // ideal case, 02/06/2014
                // Transfer bias into um, Henry Wang 03/05/2014
                // Transfer bias to nm, Henry 11/03/2015
                // The input unit for bias into ReproduceArc should be um: bias/1000
                //
                // New modification on 08/11/2016 by Henry Wang:
                // The offaxis zoneplate move towards -y axis by Rc to put the zoneplate in
                // the center of the block size.
                
                ReproduceArc(0 + 6 * l + totalslice0 * 6) = eachzoner(l) + ((float)bias / 1000) / 2;
                ReproduceArc(1 + 6 * l + totalslice0 * 6) = eachzonetheta(l);
                ReproduceArc(2 + 6 * l + totalslice0 * 6) = eachzonedr(l) - (float)bias / 1000;
                ReproduceArc(3 + 6 * l + totalslice0 * 6) = eachzonedtheta(l);
                ReproduceArc(4 + 6 * l + totalslice0 * 6) = eachzonexc(l);
                ReproduceArc(5 + 6 * l + totalslice0 * 6) = eachzoneyc(l) - Rc;
                //Leave it here in case you want to put it to the correct position,
                //and make a comparison with the parent zone plate.
                //ReproduceArc(5 + 6 * l + totalslice0 * 6) = eachzoneyc(l) ;
            }
            
            totalslice += arcNum;
            
            
        }//End of eles loop at line 1511
    }//end of for loop(zoneNum);
    
    
    if (file_format == 4 || file_format == 2) {
        
        max_dose = dosage.maxCoeff();
        min_dose = dosage.minCoeff();
        
        char  fName2[100];
        strcpy(fName2, *argv_test);
        strcat(fName2, "_ZPInfo.txt");
        string fName = fName2;
        
        FILE * fp;
        if ((fp = fopen(fName.c_str(), "wt")) == NULL)
            printf("Error");
        
        fprintf(fp, "Dose %f %f\n", max_dose, min_dose);
        
        
        for (int q = 0; q < h; q++) {
            
            double a = 1 / dosage(q);
            double b = 1 / max_dose;
            double c = 1 / min_dose;
            
            if (b != c) {
                clock = floor(((a - b) / (c - b))*(pow((long double)2, 16) - 1));
                
                //To deal with the situation that floor(trapClock) = -1;
                if (clock < 0) {
                    clock = 0;
                }
            }
            else {
                clock = (pow((long double)2, 16) - 1);
            }
            int zoneN = zoneNumber(q);
            double radius = Radius(q);
            fprintf(fp, "Zone %d %ld %f\n", zoneN, (long)clock, radius);
        }
    }
    
    // Extrac information from ReproduceArc to ReproduceArc_array for operation
    
    float* ReproduceArc_array = new float[totalslice * 6];
    totalArcs = totalslice;
    
    for (int t = 0; t < totalslice * 6; t++) {
        
        ReproduceArc_array[t] = ReproduceArc(t);
        
    }
    //Free the memory from ReproduceArc
    ReproduceArc.resize(0);
    
    // std::cout << "Here is the ZP:\n" << ReproduceArc_array << std::endl;
    std::cout << "Here is the totalslice \t= " << totalslice << std::endl;
    return ReproduceArc_array;
}


// Main function: generates sampling function to fft arc radius:\

float * HWcustomZP(double zTol, double lambda, double p, double q, double alpha, double beta, double gamma,
                   int ns, int ni, int no, int ne, double zpRad, double * orders, int numA, int &totalArcs,
                   double PC_term, double APD, double bias, int32_t zpID, int offaxis_aberration, double CRA,
                   double NA_C, int file_format, int FreeStanding, double W, double T, char ** argv_test,
                   double NoP, double IoP, int curl_on) {
    
    // secant solve parameters:\
    
    int MAXITER = 100;
    double TOLX = 0.001;
    
    int n, test0;
    double rs1, re1, rc1,
    rs2, re2, rc2,
    rs3, re3, rc3,
    dr, x0;
    
    
    int numZones = (ne - ns) / 2 + 1;
    int totalslice = 0;
    int meet_APD = 0;
    
    // For scale bar
    
    char stringBuffer[100];
    
    
    //Just for initial value
    VectorXf ReproduceArc(1);
    //VectorXd data(2);
    //int index = 0;
    
    VectorXi zoneNumber(numZones);
    VectorXd dosage(numZones);
    VectorXd Radius(numZones);
    double clock;
    int h = 0;
    double R;
    double max_dose;
    double min_dose;
    double IoP2;
    
    //IoP2 = IoP;
    
    if (NoP == IoP) {
        //If you want the last part (2nd piece for NoP = 2)
        IoP = 0;
    }
    
    
    // Fit the ZP with arc...... based on the algorism in arcfit3 & use secantsovle to
    // get the information of ZP like rs, re, rc ,dr from it.
    
    // Using 2 degree as starting point:\
    
    double totalarcNum = 360 / 2;
    
    for (int k = 0; k < numZones; k++) {
        
        //Determine the part is generated.
        //For multiple pattering purpose
        //if (NoP == 1 && IoP == 1){
        //Continue, nothing happen
        //}
        /*else if (fmod(k+1, NoP) !=  IoP){
         // Has to jump to the next part is interested!
         // Add the difference to it
         if(k == 0){
         //Lock down the first number.
         k = IoP2 - 1;
         }
         else{int Rem = fmod(k+1, NoP);
         if (Rem == 0){
         Rem = NoP;
         }
         k = k + (int) (IoP + NoP - Rem);
         }
         }
         */
        if (fmod(k + 1, NoP) != IoP) {
            //Bypass this situation
        }
        else {
            //the "}" is at the end of the function!
            
            //These terms need global definition
            
            int arcNum = 0;
            int zoneCt = k;
            int PC_on = 0;
            int APD_on = 0;
            int numofArcs = 0;
            
            VectorXd eachzoner(1);
            VectorXd eachzonedr(1);
            VectorXd eachzonetheta(1);
            VectorXd eachzonedtheta(1);
            VectorXd eachzonexc(1);
            VectorXd eachzoneyc(1);
            VectorXd no_ZP_region(1);
            
            //Start the iteration:\
            
            double theta = 0;
            double endp = (1 / totalarcNum) * 2 * M_PI;
            double endp0 = 0;
            double lastpiece = 0;
            double starting_point = 0;
            
            if (k % 5 == 0) {
                if (curl_on == 1) {
                    sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                            zpID, ( ( (float) (k + 1) ) / ((float)numZones) / 2));
                    system(stringBuffer);
                    //printf("Counting p0.5: %0.3f \n", (((float)(k + 1)) / ((float)numZones) / 2));
                    
                }
            }
            
            
            // Compute zone number::\
            
            n = ns + k * 2;
            
            // If it is in the region has PC & APD
            // We seperate APD region into 5 pieces based on angle and randomized its
            // Starting point and also distribute it uniformly
            // Updated on 06/22/2016
            // We have 2 different mode now. 1 is for the standard case and 2 is the opposite tone.
            if (FreeStanding == 1 || FreeStanding == 2) {
                // Free standing case
                // W = 0.6 dr, T= 6 dr;
                APD_on = 1;
                APD = (double)T / (W + T);
                x0 = sqrt(n * lambda * (p*q) / (p + q));
                rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                x0 = sqrt((n + 1) * lambda * (p*q) / (p + q));
                re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                
                rc1 = (rs1 + re1) / 2;
                double deltaR = abs(re1 - rs1) * 1000;
                double theta_min = W*deltaR / (rc1 * 1000);
                double empty_Num = (int)floor(2 * M_PI*(1 - APD) / theta_min);
                numofArcs = (int)empty_Num;
                if (numofArcs % 2 == 1) {
                    numofArcs = numofArcs - 1;
                }
                double no_ZP = 2 * M_PI*(1 - APD) / numofArcs;
                double spacing = 2 * M_PI*APD / numofArcs;
                //double testtt = ((double)(rand() % 360 + 1) / 360.0);
                
                //Set starting point between 1 to 360 degree
                //srand((unsigned) time(NULL));
                starting_point = 2 * M_PI*((double)(rand() % 360 + 1) / 360.0);
                theta = starting_point;
                no_ZP_region.conservativeResize(2 * numofArcs);
                no_ZP_region(0) = starting_point;
                
                int j;
                for (j = 1; j < (2 * numofArcs - 1); j = j + 2) {
                    no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                    no_ZP_region(j + 1) = no_ZP_region(j) + no_ZP;
                }
                
                if (j == (2 * numofArcs - 1)) {
                    no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                }
                
            }
            else {
                if (n >= ni && n <= no) {
                    
                    PC_on = 1;
                    
                    if (APD != 1) {
                        
                        APD_on = 1;
                        
                        x0 = sqrt(n * lambda * (p*q) / (p + q));
                        rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt((n + 1) * lambda * (p*q) / (p + q));
                        re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        rc1 = (rs1 + re1) / 2;
                        
                        double theta_min = 300 / (rc1 * 1000);
                        double empty_Num = (int)floor(2 * M_PI*(1 - APD) / theta_min);
                        numofArcs = (int)empty_Num;
                        if (numofArcs % 2 == 1) {
                            numofArcs = numofArcs - 1;
                        }
                        double no_ZP = 2 * M_PI*(1 - APD) / numofArcs;
                        double spacing = 2 * M_PI*APD / numofArcs;
                        
                        //Set starting point between 1 to 360 degree
                        //srand((unsigned) time(NULL));
                        starting_point = 2 * M_PI*((double)(rand() % 360 + 1) / 360.0);								theta = starting_point;
                        no_ZP_region.conservativeResize(2 * numofArcs);
                        no_ZP_region(0) = starting_point;
                        
                        int j;
                        for (j = 1; j < (2 * numofArcs - 1); j = j + 2) {
                            no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                            no_ZP_region(j + 1) = no_ZP_region(j) + no_ZP;
                        }
                        
                        if (j == (2 * numofArcs - 1)) {
                            no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                        }
                        
                        /*
                         no_ZP_region[1]= starting_point + spacing;
                         no_ZP_region[2]= starting_point + no_ZP + spacing;
                         no_ZP_region[3]= starting_point + no_ZP + 2*spacing;
                         no_ZP_region[4]= starting_point + 2*no_ZP + 2*spacing;
                         no_ZP_region[5]= starting_point + 2*no_ZP + 3*spacing;
                         no_ZP_region[6]= starting_point + 3*no_ZP + 3*spacing;
                         no_ZP_region[7]= starting_point + 3*no_ZP + 4*spacing;
                         no_ZP_region[8]= starting_point + 4*no_ZP + 4*spacing;
                         no_ZP_region[9]= starting_point + 4*no_ZP + 5*spacing;
                         */
                        
                        
                        //Normalized angle in between 2*M_PI
                        /*int i;
                         for (i = 0; i < 199; i++){
                         
                         if (no_ZP_region[i] > 2*M_PI)
                         no_ZP_region[i] = no_ZP_region[i] - 2*M_PI;
                         
                         }
                         */
                    }
                    
                }
            }
            
            double x1, y1,
            x2, y2,
            xm, ym;
            int second_time = 0;
            int i = 0;
            
            while (theta < (2 * M_PI + starting_point)) {
                
                //If we got APD, we have to control the angle that could be use to draw ARC
                
                if (APD_on == 1) {
                    
                    meet_APD = 0;
                    if (theta >= no_ZP_region(2 * numofArcs - 1)) {
                        break;
                    }
                    // Hit the no zp region
                    if (theta >= no_ZP_region(i + 1) && theta < no_ZP_region(i + 2)) {
                        theta = no_ZP_region(i + 2);
                        i = i + 2;
                    }
                }
                
                
                // Define parameters(test):\
                
                double meettol = 0;
                int test = 0;
                
                VectorXd Rarc(2);
                VectorXd xc(2);
                VectorXd yc(2);
                VectorXd thetada(2);
                VectorXd thetada0(2);
                VectorXd thetaa0(2);
                
                // Define the starting points
                
                double ns0 = n;
                double ne0 = n + 1;
                
                //Compute rs1 & re1:\
                
                x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                rs1 = secantSolve(x0, theta, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                re1 = secantSolve(x0, theta, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                
                rc1 = (rs1 + re1) / 2;
                
                x1 = 1000 * rc1* cos(theta);
                y1 = 1000 * rc1* sin(theta);
                
                while (1) {
                    
                    // Check if theta+ endp enter the no_ZP region
                    
                    if (APD_on == 1) {
                        if (theta >= no_ZP_region(i) && theta < no_ZP_region(i + 1)) {
                            if ((theta + endp) >= no_ZP_region(i + 1)) {
                                endp = no_ZP_region(i + 1) - theta;
                                meet_APD = 1;
                            }
                        }
                    }
                    
                    // Define the ending point:\
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs2 = secantSolve(x0, theta + endp, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re2 = secantSolve(x0, theta + endp, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc2 = (rs2 + re2) / 2;
                    
                    x2 = 1000 * rc2* cos(theta + endp);
                    y2 = 1000 * rc2* sin(theta + endp);
                    
                    // Define the centre point:\
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs3 = secantSolve(x0, theta + (endp / 2), orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re3 = secantSolve(x0, theta + (endp / 2), orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc3 = (rs3 + re3) / 2;
                    
                    xm = 1000 * rc3* cos(theta + (endp / 2));
                    ym = 1000 * rc3* sin(theta + (endp / 2));
                    
                    //Find the arc with its r, xc, yc, theta, and dtheta based on these three points:\
                    
                    double xsm = (x1 + xm) / 2;
                    double ysm = (y1 + ym) / 2;
                    
                    double xme = (xm + x2) / 2;
                    double yme = (ym + y2) / 2;
                    
                    double slopeSM = (ym - y1) / (xm - x1);
                    double slopeME = (y2 - ym) / (x2 - xm);
                    
                    xc(test) = ((slopeSM*xme - slopeME*xsm) + (yme - ysm)*slopeSM*slopeME) / (slopeSM - slopeME);
                    yc(test) = (xsm - xme + ysm*slopeSM - yme*slopeME) / (slopeSM - slopeME);
                    
                    Rarc(test) = sqrt((x1 - xc(test))*(x1 - xc(test)) + (y1 - yc(test))*(y1 - yc(test)));
                    
                    
                    // Define the coordinate based on xc yc:\
                    
                    double x_arc = x1 - xc(test);
                    double y_arc = y1 - yc(test);
                    
                    //theta: Define as the angle based on the xc & yc to the middle point of the arc:\
                    
                    if (x_arc >= 0 && y_arc >= 0) {
                        thetaa0(test) = asin(y_arc / Rarc(test));
                    }
                    else if (x_arc<0 && y_arc >= 0) {
                        thetaa0(test) = M_PI - asin(y_arc / Rarc(test));
                    }
                    else if (x_arc <= 0 && y_arc<0) {
                        thetaa0(test) = M_PI + asin(abs(y_arc / Rarc(test)));
                    }
                    else {
                        thetaa0(test) = 2 * M_PI - asin(abs(y_arc / Rarc(test)));
                    }
                    
                    //dtheta:\
                    
                    double x_arca = x2 - xc(test);
                    double y_arca = y2 - yc(test);
                    
                    if (x_arca >= 0 && y_arca >= 0) {
                        thetada0(test) = asin(y_arca / Rarc(test));
                    }
                    else if (x_arca<0 && y_arca >= 0) {
                        thetada0(test) = M_PI - asin(y_arca / Rarc(test));
                    }
                    else if (x_arca <= 0 && y_arca<0) {
                        thetada0(test) = M_PI + asin(abs(y_arca / Rarc(test)));
                    }
                    else {
                        thetada0(test) = 2 * M_PI - asin(abs(y_arca / Rarc(test)));
                    }
                    //make sure thetaa0 is always bigger than theta0:\
                    
                    
                    if (thetada0(test) > thetaa0(test)) {
                        thetada(test) = thetada0(test) - thetaa0(test);
                    }
                    else {
                        thetada(test) = thetada0(test) + (2 * M_PI - thetaa0(test));
                    }
                    
                    // Reproduce all the points in the arc and by Reproduce, ready for comparison and iteration!:\
                    
                    MatrixXd coord1(2, 11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        double th_a = thetada(test)*k / 10;
                        coord1(0, k) = Rarc(test)*cos(thetaa0(test) + th_a) + xc(test);
                        coord1(1, k) = Rarc(test)*sin(thetaa0(test) + th_a) + yc(test);
                        
                    }
                    
                    MatrixXd coord2(2, 11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        double th_r = endp*k / 10;
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        double rsr = secantSolve(x0, theta + th_r, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        double rer = secantSolve(x0, theta + th_r, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        double rcr = (rsr + rer) / 2;
                        
                        coord2(0, k) = 1000 * rcr* cos(theta + th_r);
                        coord2(1, k) = 1000 * rcr* sin(theta + th_r);
                        
                    }
                    
                    // Based on coord1 & coord2, comparing the difference with tolerance to start the iteration:\
                    
                    VectorXd output(11);
                    
                    for (int k = 0; k <11; k++) {
                        
                        output(k) = sqrt((coord1(0, k) - coord2(0, k))*(coord1(0, k) - coord2(0, k)) + (coord1(1, k) - coord2(1, k))*(coord1(1, k) - coord2(1, k)));
                        
                    }
                    
                    double result = output.maxCoeff();
                    double dth = endp / 10;
                    double rs, re;
                    double drM = 0;
                    
                    
                    for (int j = 0; j < 11; j++) {
                        
                        double thetaj = theta + j*dth;
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        rs = secantSolve(x0, thetaj, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        re = secantSolve(x0, thetaj, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        drM += (re - rs);
                    }
                    
                    double result2 = 1000 * drM / 11;
                    
                    /* 	double condition = diff(coord1, coord2) / drforArc( n, theta, endp, lambda,
                     focus, mag, zpRad, orders, numA); */
                    
                    double condition = result / result2;
                    
                    if (condition > zTol) {
                        
                        if (meettol == 1) {
                            
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            theta = theta + endp0;
                            test0 = (test == 0);
                            
                            if (theta + endp0 > (2 * M_PI + starting_point)) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = endp0;
                                lastpiece = 0;
                            }
                            
                            break;
                        }
                        
                        else {
                            endp = 0.99*endp;
                            lastpiece = 0;
                            test = (test == 0);
                        }
                    }
                    else {
                        if (lastpiece == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            test0 = test;
                            theta = theta + endp;
                            break;
                        }
                        else if (meet_APD == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            test0 = test;
                            endp0 = endp;
                            theta = theta + endp;
                            if (theta + endp0 >= (2 * M_PI + starting_point) - 0.0001) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = endp0;
                                lastpiece = 0;
                            }
                            break;
                        }
                        else {
                            endp0 = endp;
                            meettol = 1;
                            test = (test == 0);
                            
                            if (theta + 1.01*endp > (2 * M_PI + starting_point)) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = 1.01*endp;
                            }
                        }
                    }
                }//end of while 1:\
                
                
                // Switch back to um 03042014 by Henry Wang
                R = (Rarc(test0) - dr / 2) / 1000;
                eachzoner(arcNum) = (Rarc(test0) - dr / 2) / 1000;
                eachzonedr(arcNum) = dr / 1000;
                eachzonetheta(arcNum) = thetaa0(test0) + (thetada(test0) / 2);
                eachzonedtheta(arcNum) = thetada(test0);
                eachzonexc(arcNum) = xc(test0) / 1000;
                eachzoneyc(arcNum) = yc(test0) / 1000;
                
                arcNum++;
                
                //x1.conservativeResize(1+arcNum);
                //y1.conservativeResize(1+arcNum);
                //x2.conservativeResize(1+arcNum);
                //y2.conservativeResize(1+arcNum);
                //xm.conservativeResize(1+arcNum);
                //ym.conservativeResize(1+arcNum);
                
                eachzoner.conservativeResize(1 + arcNum);
                eachzonedr.conservativeResize(1 + arcNum);
                eachzonetheta.conservativeResize(1 + arcNum);
                eachzonedtheta.conservativeResize(1 + arcNum);
                eachzonexc.conservativeResize(1 + arcNum);
                eachzoneyc.conservativeResize(1 + arcNum);
                
                
                
            }//End of while theta < 2pi;:\
            
            printf("Finished zone with %d arcs.\n", arcNum);
            
            if (file_format == 4 || file_format == 2) {
                
                zoneNumber(h) = n;
                dosage(h) = (dr / 1000) / (dr / 1000 - bias / 1000);
                Radius(h) = R;
                h = h + 1;
                //zoneNumber.conservativeResize(h+1);
                //dosage.conservativeResize(h+1);
            }
            
            // Save all the information for each zones
            // We need to save 6 different parameters (rows), and all the information
            // from each zone is saved in columns:\
            
            
            int totalslice0 = totalslice;
            ReproduceArc.conservativeResize((totalslice + arcNum) * 6);
            
            for (int l = 0; l < arcNum; l++) {
                
                // Need to save the data in the form of "{R, theta, dR, dtheta, Xc, Yc }":\
                // Due to fabrication, the actual zone width should be 10nm narrower than the real situation
                // 2014/02/06
                // Transfer bias into um
                
                ReproduceArc(0 + 6 * l + totalslice0 * 6) = eachzoner(l) + (bias / 1000) / 2;
                ReproduceArc(1 + 6 * l + totalslice0 * 6) = eachzonetheta(l);
                ReproduceArc(2 + 6 * l + totalslice0 * 6) = eachzonedr(l) - bias / 1000;
                ReproduceArc(3 + 6 * l + totalslice0 * 6) = eachzonedtheta(l);
                ReproduceArc(4 + 6 * l + totalslice0 * 6) = eachzonexc(l);
                ReproduceArc(5 + 6 * l + totalslice0 * 6) = eachzoneyc(l);
            }
            
            totalslice += arcNum;
            
        }//End of the else loop for multiple patterning at line 2940
    }//end of for loop(zoneNum);
    
    
    if (file_format == 4 || file_format == 2) {
        
        max_dose = dosage.maxCoeff();
        min_dose = dosage.minCoeff();
        
        char fName2[100];
        strcpy(fName2, *argv_test);
        strcat(fName2, "_ZPInfo.txt");
        string fName = fName2;
        
        FILE * fp;
        if ((fp = fopen(fName.c_str(), "wt")) == NULL)
            printf("Error");
        
        fprintf(fp, "Dose %f %f\n", max_dose, min_dose);
        
        for (int q = 0; q < h; q++) {
            
            double a = 1 / dosage(q);
            double b = 1 / max_dose;
            double c = 1 / min_dose;
            
            
            if (b != c) {
                clock = floor(((a - b) / (c - b))*(pow((long double)2, 16) - 1));
                
                //To deal with the situation that floor(trapClock) = -1;
                if (clock < 0) {
                    clock = 0;
                }
            }
            else {
                clock = (pow((long double)2, 16) - 1);
            }
            
            
            //clock = floor(((a - b)/(c -b))*(pow((long double)2, 16)- 1));
            int zoneN = zoneNumber(q);
            double radius = Radius(q);
            fprintf(fp, "Zone %d %ld %f\n", zoneN, (long)clock, radius);
        }
    }
    
    // Extrac information from ReproduceArc to ReproduceArc_array for operation
    
    float* ReproduceArc_array = new  float[totalslice * 6];
    totalArcs = totalslice;
    
    for (int t = 0; t < totalslice * 6; t++) {
        
        ReproduceArc_array[t] = ReproduceArc(t);
        
    }
    
    //Free the memory from ReproduceArc
    ReproduceArc.resize(0);
    
    // std::cout << "Here is the ZP:\n" << ReproduceArc_array << std::endl;
    std::cout << "Here is the totalslice \t= " << totalslice << std::endl;
    return ReproduceArc_array;
}



// New main function to deal with off-axis situation.
// Build at 01/20/2017 to handle the case that we need the offaxis zoneplate to
// cross the center of its parent zonepalte
// CRA != 0
// offaxis_aberration != 0 means that the center of the aberration has been move to the center
// of the off-axis zoneplate.
// Need to add 1 more function to test if it is within the zoneplate boundary.
// We cannot adjust the angle to the proper region this time because we cannot pre-
// caculate it like the apodization or the old off - axis zoneplate condition, then move the
// angle to the proper area.


float * HWcustomZP_offaxis(double zTol, double lambda, double p, double q, double alpha, double beta, double gamma,
                           int ns, int ni, int no, int ne, double zpRad, double * orders, int numA, int &totalArcs,
                           double PC_term, double APD, double bias, int32_t zpID, int offaxis_aberration, double CRA,
                           double NA_C, int file_format, int FreeStanding, double W, double T, char ** argv_test,
                           double NoP, double IoP, int curl_on) {
    
    // secant solve parameters:\
    
    int MAXITER = 20;
    double TOLX = 0.001;
    
    int n, test0;
    double rs1, re1, rc1,
    rs2, re2, rc2,
    rs3, re3, rc3,
    dr, x0;
    
    
    int numZones = (ne - ns) / 2 + 1;
    int totalslice = 0;
    int meet_APD = 0;
    int offaxis_condition = 0;
    
    // For scale bar
    
    char stringBuffer[100];
    
    
    //Just for initial value
    VectorXf ReproduceArc(1);
    //VectorXd data(2);
    //int index = 0;
    
    VectorXi zoneNumber(1);
    VectorXd dosage(1);
    VectorXd Radius(1);
    double clock;
    int h = 0;
    double R;
    double max_dose;
    double min_dose;
    double IoP2;
    
    //IoP2 = IoP;
    
    if (NoP == IoP) {
        //If you want the last part (2nd piece for NoP = 2)
        IoP = 0;
    }
    
    
    // Fit the ZP with arc...... based on the algorism in arcfit3 & use secantsovle to
    // get the information of ZP like rs, re, rc ,dr from it.
    
    // Using 2 degree as starting point:\
    // For off axis only, much finer grid to get rid of off-axis issue.
    double totalarcNum = 3600;
    
    for (int k = 0; k < numZones; k++) {
        
        if (fmod(k + 1, NoP) != IoP) {
            //Bypass this situation
        }
        else {
            //the "}" is at the end of the function!
            
            //These terms need global definition
            
            int arcNum = 0;
            int zoneCt = k;
            int PC_on = 0;
            int APD_on = 0;
            int numofArcs = 0;
            
            VectorXd eachzoner(1);
            VectorXd eachzonedr(1);
            VectorXd eachzonetheta(1);
            VectorXd eachzonedtheta(1);
            VectorXd eachzonexc(1);
            VectorXd eachzoneyc(1);
            VectorXd no_ZP_region(1);
            
            //Start the iteration:\
            
            double theta = 0;
            double endp = (1 / totalarcNum) * 2 * M_PI;
            double endp0 = 0;
            double endp_out = 0;
            double lastpiece = 0;
            double starting_point = 0;
            
            //The free standing buttress angle
            double theta_min = 0;
            
            if (k % 5 == 0) {
                if (curl_on == 1) {
                    sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                            zpID, (((float)(k + 1)) / ((float)numZones) / 2));
                    system(stringBuffer);
                    //printf("Counting p0.5: %0.3f \n", (((float)(k + 1)) / ((float)numZones) / 2));
                    
                }
            }
            
            
            // Compute zone number::\
            
            n = ns + k * 2;
            
            //If it is in the region has PC & APD
            // We seperate APD region into 5 pieces based on angle and randomized its
            // Starting point and also distribute it uniformly
            // Updated on 06/22/2016
            // We have 2 different mode now. 1 is for the standard case and 2 is the opposite tone.
            if (FreeStanding == 1 || FreeStanding == 2) {
                // Free standing case
                // W = 0.6 dr, T= 6 dr;
                APD_on = 1;
                APD = (double)T / (W + T);
                x0 = sqrt(n * lambda * (p*q) / (p + q));
                rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                x0 = sqrt((n + 1) * lambda * (p*q) / (p + q));
                re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                
                rc1 = (rs1 + re1) / 2;
                double deltaR = abs(re1 - rs1) * 1000;
                theta_min = W*deltaR / (rc1 * 1000);
                double empty_Num = (int)floor(2 * M_PI*(1 - APD) / theta_min);
                numofArcs = (int)empty_Num;
                if (numofArcs % 2 == 1 && numofArcs != 1) {
                    numofArcs = numofArcs - 1;
                }
                double no_ZP = 2 * M_PI*(1 - APD) / numofArcs;
                double spacing = 2 * M_PI*APD / numofArcs;
                //double testtt = ((double)(rand() % 360 + 1) / 360.0);
                
                //Set starting point between 1 to 360 degree
                //srand((unsigned) time(NULL));
                starting_point = 2 * M_PI*((double)(rand() % 360 + 1) / 360.0);
                theta = starting_point;
                no_ZP_region.conservativeResize(2 * numofArcs);
                no_ZP_region(0) = starting_point;
                
                int j;
                for (j = 1; j < (2 * numofArcs - 1); j = j + 2) {
                    no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                    no_ZP_region(j + 1) = no_ZP_region(j) + no_ZP;
                }
                
                if (j == (2 * numofArcs - 1)) {
                    no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                }
                
            }
            else {
                if (n >= ni && n <= no) {
                    
                    PC_on = 1;
                    
                    if (APD != 1) {
                        
                        APD_on = 1;
                        
                        x0 = sqrt(n * lambda * (p*q) / (p + q));
                        rs1 = secantSolve(x0, theta, orders, numA, n, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt((n + 1) * lambda * (p*q) / (p + q));
                        re1 = secantSolve(x0, theta, orders, numA, n + 1, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        rc1 = (rs1 + re1) / 2;
                        
                        theta_min = 300 / (rc1 * 1000);
                        double empty_Num = (int)floor(2 * M_PI*(1 - APD) / theta_min);
                        numofArcs = (int)empty_Num;
                        if (numofArcs % 2 == 1) {
                            numofArcs = numofArcs - 1;
                        }
                        double no_ZP = 2 * M_PI*(1 - APD) / numofArcs;
                        double spacing = 2 * M_PI*APD / numofArcs;
                        
                        //Set starting point between 1 to 360 degree
                        //srand((unsigned) time(NULL));
                        starting_point = 2 * M_PI*((double)(rand() % 360 + 1) / 360.0);								theta = starting_point;
                        no_ZP_region.conservativeResize(2 * numofArcs);
                        no_ZP_region(0) = starting_point;
                        
                        int j;
                        for (j = 1; j < (2 * numofArcs - 1); j = j + 2) {
                            no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                            no_ZP_region(j + 1) = no_ZP_region(j) + no_ZP;
                        }
                        
                        if (j == (2 * numofArcs - 1)) {
                            no_ZP_region(j) = no_ZP_region(j - 1) + spacing;
                        }
                    }
                    
                }
            }
            
            double x1, y1,
            x2, y2,
            xm, ym;
            int second_time = 0;
            int i = 0;
            //The parameter to determine if arc is crossing the boundary of off-axis zoneplate
            int cross = 0;
            int arc_FirstTime = 0;
            
            while (theta < (2 * M_PI + starting_point)) {
                
                
                // For off-axis case, we check if it is within the boundary to determine whether
                // we should save this piece of arc or not, or we have to adjust its starting or ending angle!
                // If both points are within the boundary => this arc is in the off-axis zoneplate
                // If Head is in but not the tail, save the previous one that when tail is still in
                // This piece has to at least larger the size of buttress to avoid issues like generate
                // arc that is totally not part of the target off-axis zoneplate
                // If Head is not in but the tail is, use the current one as the new head for the calculation
                // and don't save anything
                // You have to check if HeadIsIn or TailIsIn every time the theta and endp is changing!
                // If the status is changed you always have to save the old one since that would be the
                // new starting point or boundary you have to use for your arc!
                
                
                //If we got APD, we have to control the angle that could be use to draw ARC
                
                if (APD_on == 1) {
                    
                    meet_APD = 0;
                    if (theta >= no_ZP_region(2 * numofArcs - 1)) {
                        break;
                    }
                    // Hit the no zp region
                    if (theta >= no_ZP_region(i + 1) && theta < no_ZP_region(i + 2)) {
                        theta = no_ZP_region(i + 2);
                        i = i + 2;
                    }
                }
                
                
                // Define parameters(test):\
                
                double meettol = 0;
                int test = 0;
                int TailIdx = 0;
                
                VectorXd Rarc(2);
                VectorXd xc(2);
                VectorXd yc(2);
                VectorXd thetada(2);
                VectorXd thetada0(2);
                VectorXd thetaa0(2);
                
                //Only has 1 head but the tail is varying
                //so you need to save the old and new one
                std::vector<int> HeadIsIn(1,0);
                std::vector<int> TailIsIn(2, 0);
                
                // Define the starting points
                
                double ns0 = n;
                double ne0 = n + 1;
                
                //Compute rs1 & re1:\
                
                x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                rs1 = secantSolve(x0, theta, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                re1 = secantSolve(x0, theta, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                
                rc1 = (rs1 + re1) / 2;
                
                x1 = 1000 * rc1* cos(theta);
                y1 = 1000 * rc1* sin(theta);
                
                //Check the head of the arc is in the off-axis zoneplate
                //Head stays the same for each arc calculation
                //Only the tail would vary and need to save the old and new one
                HeadIsIn[0] = IsOffAxis(rc1, theta, p, q, CRA, NA_C);
                
                
                while (1) {
                    
                    // Check if theta+ endp enter the no_ZP region
                    
                    if (APD_on == 1) {
                        if (theta >= no_ZP_region(i) && theta < no_ZP_region(i + 1)) {
                            if ((theta + endp) >= no_ZP_region(i + 1)) {
                                endp = no_ZP_region(i + 1) - theta;
                                meet_APD = 1;
                            }
                        }
                    }
                    
                    // Define the ending point:\
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs2 = secantSolve(x0, theta + endp, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re2 = secantSolve(x0, theta + endp, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc2 = (rs2 + re2) / 2;
                    
                    //Check the tail of the arc is in the off-axis zoneplate
                    //Head stays the same for each arc calculation
                    //Only the tail would vary and need to save the old and new one
                    //First time for this arc to estimate if the Tail is in
                    TailIsIn[TailIdx] = IsOffAxis(rc2, theta + endp, p, q, CRA, NA_C);
                    TailIdx = !TailIdx;
                    
                    //Has to check it "after" there is at least one set of arc information!!!
                    // Has to make decision here to see you need to cut it off
                    if (arc_FirstTime != 0) {
                        if (HeadIsIn[0] == 0)
                        {
                            if (TailIsIn[(!TailIdx)] != 0) {
                                //This piece doesn't belong to the offaxis zoneplate
                                //Use the previous endp and cut the arc right here
                                //Use the average of the current endp + theta as new starting point
                                //for theta
                                dr = drforArc(n, theta, endp / 1.01, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                                endp_out = endp / 1.01;
                                theta = theta + endp;
                                endp = (1 / totalarcNum) * 2 * M_PI;
                                cross = 1;
                                //You have to use the old set of arc information for the ouput
                                //Since you check arc_FirstTime already, you know there is for sure
                                //Have a set of data that is ready for output
                                test0 = !test;
                                offaxis_condition = 0;
                                break;
                            }
                        }
                        else if (HeadIsIn[0] == 1)
                        {
                            if (TailIsIn[(!TailIdx)] != 1) {
                                //This piece belong to the offaxis zonepalte
                                //Use the previous endp and cut the arc right here
                                //Use the average of the current and old endp + theta as new starting point
                                //for theta
                                dr = drforArc(n, theta, endp / 1.01, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                                theta = theta + endp;
                                endp_out = endp / 1.01;
                                //Has to check if this piece is too small to keep in the target off-axis zoneplate
                                //which might generate unwanted pattern that is not the target off-axis zoneplate at all.
                                endp = (1 / totalarcNum) * 2 * M_PI;
                                cross = 1;
                                //You have to use the old set of arc information for the ouput
                                //Since you check arc_FirstTime already, you know there is for sure
                                //Have a set of data that is ready for output
                                test0 = !test;
                                
                                offaxis_condition = 1;
                                
                                
                                break;
                            }
                        }
                    }
                    
                    x2 = 1000 * rc2* cos(theta + endp);
                    y2 = 1000 * rc2* sin(theta + endp);
                    
                    // Define the centre point:\
                    
                    x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                    rs3 = secantSolve(x0, theta + (endp / 2), orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                    re3 = secantSolve(x0, theta + (endp / 2), orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                    
                    rc3 = (rs3 + re3) / 2;
                    
                    xm = 1000 * rc3* cos(theta + (endp / 2));
                    ym = 1000 * rc3* sin(theta + (endp / 2));
                    
                    //Find the arc with its r, xc, yc, theta, and dtheta based on these three points:\
                    
                    double xsm = (x1 + xm) / 2;
                    double ysm = (y1 + ym) / 2;
                    
                    double xme = (xm + x2) / 2;
                    double yme = (ym + y2) / 2;
                    
                    double slopeSM = (ym - y1) / (xm - x1);
                    double slopeME = (y2 - ym) / (x2 - xm);
                    
                    xc(test) = ((slopeSM*xme - slopeME*xsm) + (yme - ysm)*slopeSM*slopeME) / (slopeSM - slopeME);
                    yc(test) = (xsm - xme + ysm*slopeSM - yme*slopeME) / (slopeSM - slopeME);
                    
                    Rarc(test) = sqrt((x1 - xc(test))*(x1 - xc(test)) + (y1 - yc(test))*(y1 - yc(test)));
                    
                    
                    // Define the coordinate based on xc yc:\
                    
                    double x_arc = x1 - xc(test);
                    double y_arc = y1 - yc(test);
                    
                    //theta: Define as the angle based on the xc & yc to the middle point of the arc:\
                    
                    if (x_arc >= 0 && y_arc >= 0) {
                        thetaa0(test) = asin(y_arc / Rarc(test));
                    }
                    else if (x_arc<0 && y_arc >= 0) {
                        thetaa0(test) = M_PI - asin(y_arc / Rarc(test));
                    }
                    else if (x_arc <= 0 && y_arc<0) {
                        thetaa0(test) = M_PI + asin(abs(y_arc / Rarc(test)));
                    }
                    else {
                        thetaa0(test) = 2 * M_PI - asin(abs(y_arc / Rarc(test)));
                    }
                    
                    //dtheta:\
                    
                    double x_arca = x2 - xc(test);
                    double y_arca = y2 - yc(test);
                    
                    if (x_arca >= 0 && y_arca >= 0) {
                        thetada0(test) = asin(y_arca / Rarc(test));
                    }
                    else if (x_arca<0 && y_arca >= 0) {
                        thetada0(test) = M_PI - asin(y_arca / Rarc(test));
                    }
                    else if (x_arca <= 0 && y_arca<0) {
                        thetada0(test) = M_PI + asin(abs(y_arca / Rarc(test)));
                    }
                    else {
                        thetada0(test) = 2 * M_PI - asin(abs(y_arca / Rarc(test)));
                    }
                    //make sure thetaa0 is always bigger than theta0:\
                    
                    
                    if (thetada0(test) > thetaa0(test)) {
                        thetada(test) = thetada0(test) - thetaa0(test);
                    }
                    else {
                        thetada(test) = thetada0(test) + (2 * M_PI - thetaa0(test));
                    }
                    
                    // Reproduce all the points in the arc and by Reproduce, ready for comparison and iteration!:\
                    
                    MatrixXd coord1(2, 11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        double th_a = thetada(test)*k / 10;
                        coord1(0, k) = Rarc(test)*cos(thetaa0(test) + th_a) + xc(test);
                        coord1(1, k) = Rarc(test)*sin(thetaa0(test) + th_a) + yc(test);
                        
                    }
                    
                    MatrixXd coord2(2, 11);
                    
                    for (int k = 0; k < 11; k++) {
                        
                        double th_r = endp*k / 10;
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        double rsr = secantSolve(x0, theta + th_r, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        double rer = secantSolve(x0, theta + th_r, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        double rcr = (rsr + rer) / 2;
                        
                        coord2(0, k) = 1000 * rcr* cos(theta + th_r);
                        coord2(1, k) = 1000 * rcr* sin(theta + th_r);
                        
                    }
                    
                    // Based on coord1 & coord2, comparing the difference with tolerance to start the iteration:\
                    
                    VectorXd output(11);
                    
                    for (int k = 0; k <11; k++) {
                        
                        output(k) = sqrt((coord1(0, k) - coord2(0, k))*(coord1(0, k) - coord2(0, k)) + (coord1(1, k) - coord2(1, k))*(coord1(1, k) - coord2(1, k)));
                        
                    }
                    
                    double result = output.maxCoeff();
                    double dth = endp / 10;
                    double rs, re;
                    double drM = 0;
                    
                    
                    for (int j = 0; j < 11; j++) {
                        
                        double thetaj = theta + j*dth;
                        x0 = sqrt(ns0 * lambda * (p*q) / (p + q));
                        rs = secantSolve(x0, thetaj, orders, numA, ns0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        x0 = sqrt(ne0 * lambda * (p*q) / (p + q));
                        re = secantSolve(x0, thetaj, orders, numA, ne0, p, lambda, q, alpha, beta, gamma, zpRad, 1, MAXITER, TOLX, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                        
                        drM += (re - rs);
                    }
                    
                    double result2 = 1000 * drM / 11;
                    
                    /* 	double condition = diff(coord1, coord2) / drforArc( n, theta, endp, lambda,
                     focus, mag, zpRad, orders, numA); */
                    
                    double condition = result / result2;
                    
                    if (condition > zTol) {
                        
                        if (meettol == 1) {
                            
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            theta = theta + endp0;
                            //test0 is the old test, in order to retrieve the previous set of data
                            test0 = !test;
                            endp_out = endp;
                            
                            //endp = (1 / totalarcNum) * 2 * M_PI
                            
                            if (theta + ((1 / totalarcNum) * 2 * M_PI) > (2 * M_PI + starting_point)) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                endp = (1 / totalarcNum) * 2 * M_PI;
                                lastpiece = 0;
                            }
                            // The current TailIdx points to the previous one, not the current tail angle!
                            break;
                        }
                        
                        else {
                            //Fix the step size in order to fix the boundary issue
                            // at large endp, 0.01 is a huge step
                            endp = 0.99* endp;
                            lastpiece = 0;
                            //Rotating the index
                            test = !test;
                        }
                    }
                    else {
                        //For checking off-axis boundary condition
                        //It has to have at least 1 set of valid arc information
                        //Before you make the decision to cut it or not due to the boundary condition.
                        arc_FirstTime = 1;
                        
                        if (lastpiece == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            //The current idx is the idx used for output
                            test0 = test;
                            theta = theta + endp;
                            endp_out = endp;
                            break;
                        }
                        else if (meet_APD == 1) {
                            dr = drforArc(n, theta, endp, lambda, p, q, alpha, beta, gamma, zpRad, orders, numA, PC_term, PC_on, offaxis_aberration, CRA, NA_C);
                            //The current idx is the idx used for output
                            test0 = test;
                            endp0 = endp;
                            endp_out = endp;
                            theta = theta + endp;
                            if (theta + ((1 / totalarcNum) * 2 * M_PI) >= (2 * M_PI + starting_point) - 0.0001) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                //endp = (1 / totalarcNum) * 2 * M_PI
                                //Always start with the smaller piece
                                endp = (1 / totalarcNum) * 2 * M_PI;
                                lastpiece = 0;
                            }
                            // The current TailIdx points to the previous one, not the current tail angle!
                            
                            break;
                        }
                        else {
                            endp0 = endp;
                            meettol = 1;
                            //rotate the idx
                            test = !test;
                            
                            if (theta + 1.01*endp > (2 * M_PI + starting_point)) {
                                endp = (2 * M_PI + starting_point) - theta;
                                lastpiece = 1;
                            }
                            else {
                                //Fix the step size in order to fix the boundary issue
                                // at large endp, 0.01 is a huge step
                                endp = endp + (1 / totalarcNum) * 2 * M_PI;
                            }
                        }
                    }
                }//end of while 1:\
                
                
                //Has to setup a condition to check the arc condition
                //cross == 0 : stay in either inside or outside the offaxis region
                //             just make simple judgement to see if we want to keep this arc
                //cross != 1 : crossing the boundary, depends on ethier tail is out or tail is in
                //             special situation, has to deal with it seperately.
                if (cross == 0) {
                    if (HeadIsIn[0] == 1) {
                        if (TailIsIn[!TailIdx] == 1)
                        {//This piece is belong to the off-axis zoneplate
                            offaxis_condition = 1;
                        }
                    }
                    else
                    {
                        //If the current head is not in there
                        if (TailIsIn[!TailIdx] == 0)
                        {//This piece is not belong to the off-axis zoneplate
                            offaxis_condition = 0;
                        }
                    }
                }
                
                
                
                
                if (offaxis_condition == 1){
                    
                    //Make sure that every piece get output is at
                    //least in the same size as the buttress
                    //Add 0.9* to further increase the limitation
                    //Trying to avoid any piece which is smaller than half the unit endp
                    //The piece has to be at least larger then the butterss (0.6*dr)
                    //The size of the arc (r_arc * dtheta) has to larger than the width of the
                    //butterss (0.6*dr)
                    if ((((Rarc(test0) - dr / 2) / 1000)* thetada(test0)) < (0.6 * dr / 1000)) {
                        //offaxis_condition = 0;
                        //just don't save the arc
                    }
                    else{
                        //Only save the piece that is belong to the offaxis zoneplate
                        // Switch back to um 03042014 by Henry Wang
                        
                        // Deal with extreme situation, only has 1 point
                        // in the off-axis zoneplate, never actually get to
                        // the formation of arc
                        if (Rarc(test0) >= 0) {
                            R = (Rarc(test0) - dr / 2) / 1000;
                            eachzoner(arcNum) = (Rarc(test0) - dr / 2) / 1000;
                            eachzonedr(arcNum) = dr / 1000;
                            eachzonetheta(arcNum) = thetaa0(test0) + (thetada(test0) / 2);
                            eachzonedtheta(arcNum) = thetada(test0);
                            eachzonexc(arcNum) = xc(test0) / 1000;
                            eachzoneyc(arcNum) = yc(test0) / 1000;
                            
                            arcNum++;
                            
                            eachzoner.conservativeResize(1 + arcNum);
                            eachzonedr.conservativeResize(1 + arcNum);
                            eachzonetheta.conservativeResize(1 + arcNum);
                            eachzonedtheta.conservativeResize(1 + arcNum);
                            eachzonexc.conservativeResize(1 + arcNum);
                            eachzoneyc.conservativeResize(1 + arcNum);
                        }
                    }
                }
                
                //Reset the condition for the next arc
                offaxis_condition = 0;
                cross = 0;
                arc_FirstTime = 0;
                
            }//End of while theta < 2pi;:\
            
            printf("Finished zone with %d arcs.\n", arcNum);
            
            //We only keep track of dosage for those actually generate the
            //arc for offaxis zoneplate
            if (arcNum != 0) {
                
                if (file_format == 4 || file_format == 2) {
                    
                    if (h != 0) {
                        zoneNumber.conservativeResize(h + 1);
                        dosage.conservativeResize(h + 1);
                        Radius.conservativeResize(h + 1);
                    }
                    zoneNumber(h) = n;
                    if (n == 1239)
                    {
                        int stop = 0;
                    }
                    dosage(h) = (dr / 1000) / (dr / 1000 - bias / 1000);
                    Radius(h) = R;
                    h = h + 1;
                    
                }
            }
            
            // Save all the information for each zones
            // We need to save 6 different parameters (rows), and all the information
            // from each zone is saved in columns:\
										  
            
            int totalslice0 = totalslice;
            ReproduceArc.conservativeResize((totalslice + arcNum) * 6);
            
            for (int l = 0; l < arcNum; l++) {
                
                // Need to save the data in the form of "{R, theta, dR, dtheta, Xc, Yc }":\
                // Due to fabrication, the actual zone width should be 10nm narrower than the real situation
                // 2014/02/06
                // Transfer bias into um
                
                ReproduceArc(0 + 6 * l + totalslice0 * 6) = eachzoner(l) + (bias / 1000) / 2;
                ReproduceArc(1 + 6 * l + totalslice0 * 6) = eachzonetheta(l);
                ReproduceArc(2 + 6 * l + totalslice0 * 6) = eachzonedr(l) - bias / 1000;
                ReproduceArc(3 + 6 * l + totalslice0 * 6) = eachzonedtheta(l);
                ReproduceArc(4 + 6 * l + totalslice0 * 6) = eachzonexc(l);
                ReproduceArc(5 + 6 * l + totalslice0 * 6) = eachzoneyc(l);
            }
            
            totalslice += arcNum;
            
        }//End of the else loop for multiple patterning at line 2940
    }//end of for loop(zoneNum);
    
    
    if (file_format == 4 || file_format == 2) {
        
        max_dose = dosage.maxCoeff();
        min_dose = dosage.minCoeff();
        
        char fName2[100];
        strcpy(fName2, *argv_test);
        strcat(fName2, "_ZPInfo.txt");
        string fName = fName2;
        
        FILE * fp;
        if ((fp = fopen(fName.c_str(), "wt")) == NULL)
            printf("Error");
        
        fprintf(fp, "Dose %f %f\n", max_dose, min_dose);
        
        for (int q = 0; q < h; q++) {
            
            double a = 1 / dosage(q);
            double b = 1 / max_dose;
            double c = 1 / min_dose;
            
            
            if (b != c) {
                clock = floor(((a - b) / (c - b))*(pow((long double)2, 16) - 1));
                
                //To deal with the situation that floor(trapClock) = -1;
                if (clock < 0) {
                    clock = 0;
                }
            }
            else {
                clock = (pow((long double)2, 16) - 1);
            }
            
            
            //clock = floor(((a - b)/(c -b))*(pow((long double)2, 16)- 1));
            int zoneN = zoneNumber(q);
            double radius = Radius(q);
            fprintf(fp, "Zone %d %ld %f\n", zoneN, (int32_t) clock, radius);
        }
    }
    
    // Extrac information from ReproduceArc to ReproduceArc_array for operation
    
    float* ReproduceArc_array = new  float[totalslice * 6];
    totalArcs = totalslice;
    
    for (int t = 0; t < totalslice * 6; t++) {
        
        ReproduceArc_array[t] = ReproduceArc(t);
        
    }
    
    //Free the memory from ReproduceArc
    ReproduceArc.resize(0);
    
    // std::cout << "Here is the ZP:\n" << ReproduceArc_array << std::endl;
    std::cout << "Here is the totalslice \t= " << totalslice << std::endl;
    return ReproduceArc_array;
}

// writes to text file :

void writeToGDSCell(float* xArray, double arrayR, double arrayC, char** argv_test, int32_t zpID, int curl_on) {
    
    string dateStr = "dateString";
    string fName = *argv_test;
    string cName = "cellName";
    char    stringBuffer[100];
    
    FILE * fp;
    if ((fp = fopen(fName.c_str(), "wt")) == NULL)
        printf("Error");
    
    for (int k = 0; k < arrayR; k++) {
        if (curl_on == 1) {
            if (k % 5 == 0) {
                sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                        zpID, ( ( (float) (k+1) ) / ((float) arrayR) / 2) + 0.5);
                system(stringBuffer);
                //printf("Counting m0.5: %0.3f \n", (((float)(k + 1)) / ((float)arrayR) / 2) + 0.5);
                
            }
        }
        
        for (int m = 0; m < arrayC; m++) {
            fprintf(fp, "%0.4f", xArray[6 * k + m]);
            if (m != arrayC - 1)
                fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

void writeToGDSCell_Trap(float* xArray, int arrayR, int arrayC, char** argv_test, int32_t zpID, double first, FILE* fp, int curl_on) {
    
    //FILE * fp;
    char    stringBuffer[100];
    
    //First time then open the file
    if (first == 0) {
        //    string dateStr  = "2013-12-12 00:00:00";
        string fName = *argv_test;
        //    string cName    = "cellName";
        //FILE * fp;
        if ((fp = fopen(fName.c_str(), "wt")) == NULL)
            printf("Error");
    }
    
    //		fprintf(fp, "gds2{7\nm=%s a=%s\nlib 'um.DB' 1e-04 1e-10\n", dateStr.c_str(), dateStr.c_str());
    //fprintf(fp, "gds2{3\nm=%s a=%s\nlib 'noname' 1 1e-09\n", dateStr.c_str(), dateStr.c_str());
    //		fprintf(fp, "cell{c=%s m=%s '%s'\n", dateStr.c_str(), dateStr.c_str(), cName.c_str());
    
    for (int k = 0; k < arrayR; k++) {
        //			fprintf(fp, "b{1 xy(");
        
        for (int m = 0; m < arrayC; m++) {
            fprintf(fp, "%0.4f", xArray[arrayR*m + k]);
            if (m != arrayC - 1)
                fprintf(fp, " ");
        }
        fprintf(fp, "\n");
        //        fprintf(fp, ")}\n");
    }
    
    
    //Close the file
    if (first == 1) {
        //		fprintf(fp, "}\n}");
        fprintf(fp, "\n");
        //	fclose(fp);
    }
    
}

unsigned long decode32(unsigned char * int32) {
    
    /*for (int k = 0; k < 4; k++){
     printf("%d\n", (int) int32[k]);
     }*/
    unsigned long out = 0;
    
    out += ((unsigned long)int32[0]) << 0;
    out += ((unsigned long)int32[1]) << 8;
    out += ((unsigned long)int32[2]) << 16;
    out += ((unsigned long)int32[3]) << 24;
    return out;
}

void encode32(long coord, unsigned char * cPart) {
    
    cPart[0] = (coord >> 24) & 255;
    cPart[1] = (coord >> 16) & 255;
    cPart[2] = (coord >> 8) & 255;
    cPart[3] = (coord)& 255;
}

void encode16(long val, unsigned char * cPart) {
    cPart[0] = (val >> 8) & 255;
    cPart[1] = (val)& 255;
}


// Gateway function
void mexFunction(float* xArray, int arrayR, int arrayC, char** argv_test, int32_t zpID,
                 double first, FILE* op, int name_size, int curl_on)
{
    //first is used to check if it is the vary first piece or the last one
    //To create or close the file
    //first = 0: 1st time enter the function
    //first = 1: last time enter the function, close the file
    // arrayR = 0 since we deal with one trap at a time
    
    unsigned char gdsPost[8];
    unsigned char polyPre[16];
    unsigned char polyPost[4];
    unsigned char polyForm[2];
    //FILE *op;
    
    int* gds_Name_preamble = NULL;
    
    char *buffer;
    unsigned char sizeInfo[4];
    
    //Write gds file
    // Read number of polygons:
    //long nPoly = arrayR;
    // Define ambles:
    
    int gdspostamble[8] = { 0, 4, 7, 0, 0, 4, 4, 0 };
    int polypreamble[16] = { 0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0 };
    int polypostamble[4] = { 0, 4, 17, 0 };
    
    int polyBlockFormat[2] = { 16, 3 };
    
    for (int k = 0; k < 8; k++)
        //gdsPost[k] = (unsigned char)gdspostamble[k];
        *(gdsPost + k) = (unsigned char)gdspostamble[k];
    for (int k = 0; k < 16; k++)
        //polyPre[k] = (unsigned char)polypreamble[k];
        *(polyPre + k) = (unsigned char)polypreamble[k];
    for (int k = 0; k < 4; k++)
        *(polyPost + k) = (unsigned char)polypostamble[k];
    for (int k = 0; k < 2; k++)
        *(polyForm + k) = (unsigned char)polyBlockFormat[k];
    
    
    //unsigned char buffer2[4];
    unsigned char polyFormNBuffer[2]; // stores number of bytes
    
    //unsigned long numVertices;
    //long numBytes;
    
    for (int k = 0; k < arrayR; k++) {
        
        //numVertices = 4;
        //numBytes =  (4 * 2 * 4);
        
        
        
        // write polypreamble
        fwrite(polyPre, sizeof(char), 16, op);
        
        // write polyform:
        encode16(44, polyFormNBuffer);
        //double a = decode32(polyFormNBuffer);
        
        //printf("numBytes: [%d]\n", (int)numBytes);
        //printf("polyFormNBuffer: [%d, %d]\n", (int)polyFormNBuffer[0], (int)polyFormNBuffer[1]);
        
        fwrite(polyFormNBuffer, sizeof(char), 2, op);
        fwrite(polyForm, sizeof(char), 2, op);
        
        // read and write block
        unsigned char dataBuffer2[40];
        unsigned char dataBuffer[4];
        
        
        for (int m = 0; m < arrayC; m++) {
            int32_t num = xArray[arrayR*m + k] * 10000;
            encode32(num, dataBuffer);
            for (int q = 0; q < 4; q++) {
                dataBuffer2[m * 4 + q] = dataBuffer[q];
            }
        }
        
        //Rewrite the 1st verticies
        for (int m = 0; m < 2; m++) {
            int32_t num = xArray[arrayR*m + k] * 10000;
            encode32(num, dataBuffer);
            for (int q = 0; q< 4; q++) {
                dataBuffer2[(m + 8) * 4 + q] = dataBuffer[q];
            }
        }
        
        
        
        fwrite(dataBuffer2, 1, 40, op);
        
        //thisPolyBlock += polypostamble;
        fwrite(polyPost, sizeof(char), 4, op);
    }//end of arrayR
    
    if (first == 1) {
        // Write GDS postamble:
        fwrite(gdsPost, sizeof(char), 8, op);
        //fclose(op);
    }
    
    //free(xArray);
    //free(gds_Name_preamble);
    //free(polyPre);
    //free(polyForm);
    //free(polyPost);
    //free(gdsPost);
    
}


// Read 8 vertices of each polygons, and then cut it into 2 triangle and 1 rectangle.
void CAD_tool(float* xArray, int numTrap, int arrayC, char** argv_test, int32_t zpID, int block_size, int wrv_split,
              int block_index, int32_t & totalPoly, double maxDose, double minDose, double first, double* first_WRV_Split, FILE* fp, FILE** fp_wrv) {
    //For trap, the format will be x1 y1 x3 y3 x2 x4
    // x1 y1 is the point at bottom left, x3 y3 is at the top rihgt.
    // x2 is at the bottom right and if x2 = x1 then it is triangle
    // x4 is at the top left and if x4 = x3 then it is triangle
    // Include a -90 degree rotation matrix to transform the zp to the +x-axis
    // Transform first, then smart cut, and then add the displacement
    // x' = y
    // y' = -x
    // Write to the file directly to save the memory allocation problem
    //
    //06/24/2016
    //Adding the ability to handle 9 blocks together if the wrv_split = 1
    //Choose the right fp before get into this one
    //New index to indentify the block index
    
    
    //Take the coordiinates into cacluation, then sort it out
    VectorXf vertices_coord(8);
    VectorXf vertices_coord3(8);
    VectorXd SmartCut;
    int32_t vertices_coord2[8];
    int32_t vertices_test[8];
    int32_t vertices_dose;
    VectorXf vertices_ycoord(4);
    int32_t vertices_ycoord_test[4];
    //Rectangle and triangle
    int32_t Trap_Rec[6];
    int32_t Trap_Tri1[6];
    int32_t Trap_Tri2[6];
    int32_t Smart_Cut[6];
    int32_t Smart_Cut2[7];
    VectorXf one_set2(7);
    totalPoly = 0;
    int32_t AfterCut[15];
    int32_t Trap[8];
    int32_t Tri[6];
    int32_t vertices_coord2_y[4];
    int32_t extra_x, extra_y;
    int test_situation;
    //Unit nm
    double pixel_size = 0.5;
    FILE* outputFile = NULL;
    std::string fName;
    if (wrv_split == 1) {
        outputFile = fp_wrv[block_index - 1];
    }
    
    //string fName = *argv_test;
    //FILE * fp;
    
    if (wrv_split == 1) {
        //If it has wrv split
        //You need to check the individual first to see if you need to do it again.
        //If we split the wrv file, we have to speficy now the block size is only 80% of its original
        //Thus the vepdef has to change to 80% of its original value.
        
        if (first_WRV_Split[block_index - 1] == 0) {
            fprintf(outputFile, "patdef 500 %d %d %f %f 0 0\n", block_size / 10, block_size / 10, maxDose, minDose);
            fprintf(outputFile, "vepdef 20 %d %d\n", (int)(block_size / 10 * 0.8), (int)(block_size / 10 * 0.8));
        }
    }
    else {
        // Only 1 first
        if (first == 0) {
            //if ((fp = fopen(fName.c_str(), "wt")) == NULL)
            //	printf("Error");
            
            fprintf(fp, "patdef 500 %d %d %f %f 0 0\n", block_size / 10, block_size / 10, maxDose, minDose);
            fprintf(fp, "vepdef 20 %d %d\n", block_size / 10, block_size / 10);
        }
    }
    // The offset of the center has to be modified based on resolution unit
    int new_center = block_size / 2 / 2;
    int num_1st, num_2nd, num_3rd, num_4th, index_1st, index_2nd, index_3rd, index_4th;
    int pass = 0;
    int piece;
    
    // arrayR represent numberTrap, now for smart cut,
    // Each trapezoid split into 2 triangle and 1 rectangle
    // Total vertices= (2*3+1*4)*2 = 20 vertices 
    // Need conservative resize
    //VectorXd CAD_smartcut[1];
    
    
    for (int k = 0; k < numTrap; k++) {
        
        //printf("%lld out of total trap %lld \n", k+1, numTrap+1);
        //printf("Total polygon %ld \n", totalPoly);		
        
        //Turn the coordinates from um to A.
        for (int m = 0; m < arrayC; m++) {
            vertices_coord(m) = xArray[numTrap*m + k] * 10000;
        }
        vertices_dose = floor(xArray[numTrap * 8 + k]);
        
        //Do the rotation first: -90 degree
        vertices_coord3(0) = vertices_coord(1) + new_center;
        vertices_coord3(1) = -1 * vertices_coord(0) + new_center;
        vertices_coord3(2) = vertices_coord(3) + new_center;
        vertices_coord3(3) = -1 * vertices_coord(2) + new_center;
        vertices_coord3(4) = vertices_coord(5) + new_center;
        vertices_coord3(5) = -1 * vertices_coord(4) + new_center;
        vertices_coord3(6) = vertices_coord(7) + new_center;
        vertices_coord3(7) = -1 * vertices_coord(6) + new_center;
        
        
        if (wrv_split == 1) {
            // Based on the wrv_split and block number
            // Move to the correct position in the corresponding block with new center.
            if (block_index == 1) {
                // x' = x + 400, y' = y - 400
                vertices_coord3(0) = vertices_coord3(0) + 4000000;
                vertices_coord3(1) = vertices_coord3(1) - 4000000;
                vertices_coord3(2) = vertices_coord3(2) + 4000000;
                vertices_coord3(3) = vertices_coord3(3) - 4000000;
                vertices_coord3(4) = vertices_coord3(4) + 4000000;
                vertices_coord3(5) = vertices_coord3(5) - 4000000;
                vertices_coord3(6) = vertices_coord3(6) + 4000000;
                vertices_coord3(7) = vertices_coord3(7) - 4000000;
            }
            else if (block_index == 2) {
                // x' = x, y' = y - 400
                vertices_coord3(0) = vertices_coord3(0);
                vertices_coord3(1) = vertices_coord3(1) - 4000000;
                vertices_coord3(2) = vertices_coord3(2);
                vertices_coord3(3) = vertices_coord3(3) - 4000000;
                vertices_coord3(4) = vertices_coord3(4);
                vertices_coord3(5) = vertices_coord3(5) - 4000000;
                vertices_coord3(6) = vertices_coord3(6);
                vertices_coord3(7) = vertices_coord3(7) - 4000000;
                
            }
            else if (block_index == 3) {
                // x' = x - 400, y' = y - 400
                vertices_coord3(0) = vertices_coord3(0) - 4000000;
                vertices_coord3(1) = vertices_coord3(1) - 4000000;
                vertices_coord3(2) = vertices_coord3(2) - 4000000;
                vertices_coord3(3) = vertices_coord3(3) - 4000000;
                vertices_coord3(4) = vertices_coord3(4) - 4000000;
                vertices_coord3(5) = vertices_coord3(5) - 4000000;
                vertices_coord3(6) = vertices_coord3(6) - 4000000;
                vertices_coord3(7) = vertices_coord3(7) - 4000000;
                
            }
            else if (block_index == 4) {
                // Move to the new origin, the coordinate has to do the 
                // following
                // x' = x - (-400), y' = y
                vertices_coord3(0) = vertices_coord3(0) + 4000000;
                vertices_coord3(1) = vertices_coord3(1);
                vertices_coord3(2) = vertices_coord3(2) + 4000000;
                vertices_coord3(3) = vertices_coord3(3);
                vertices_coord3(4) = vertices_coord3(4) + 4000000;
                vertices_coord3(5) = vertices_coord3(5);
                vertices_coord3(6) = vertices_coord3(6) + 4000000;
                vertices_coord3(7) = vertices_coord3(7);
                
            }
            else if (block_index == 5) {
                //The center (original) block
                //Nothing need to be changed
                //x' = x, y' = y
            }
            else if (block_index == 6) {
                // x' = x - 400, y' = y
                vertices_coord3(0) = vertices_coord3(0) - 4000000;
                vertices_coord3(1) = vertices_coord3(1);
                vertices_coord3(2) = vertices_coord3(2) - 4000000;
                vertices_coord3(3) = vertices_coord3(3);
                vertices_coord3(4) = vertices_coord3(4) - 4000000;
                vertices_coord3(5) = vertices_coord3(5);
                vertices_coord3(6) = vertices_coord3(6) - 4000000;
                vertices_coord3(7) = vertices_coord3(7);
                
            }
            else if (block_index == 7) {
                // x' = x + 400, y' = y + 400
                vertices_coord3(0) = vertices_coord3(0) + 4000000;
                vertices_coord3(1) = vertices_coord3(1) + 4000000;
                vertices_coord3(2) = vertices_coord3(2) + 4000000;
                vertices_coord3(3) = vertices_coord3(3) + 4000000;
                vertices_coord3(4) = vertices_coord3(4) + 4000000;
                vertices_coord3(5) = vertices_coord3(5) + 4000000;
                vertices_coord3(6) = vertices_coord3(6) + 4000000;
                vertices_coord3(7) = vertices_coord3(7) + 4000000;
            }
            else if (block_index == 8) {
                // x' = x, y' = y + 400
                vertices_coord3(0) = vertices_coord3(0);
                vertices_coord3(1) = vertices_coord3(1) + 4000000;
                vertices_coord3(2) = vertices_coord3(2);
                vertices_coord3(3) = vertices_coord3(3) + 4000000;
                vertices_coord3(4) = vertices_coord3(4);
                vertices_coord3(5) = vertices_coord3(5) + 4000000;
                vertices_coord3(6) = vertices_coord3(6);
                vertices_coord3(7) = vertices_coord3(7) + 4000000;
            }
            else if (block_index == 9) {
                // x' = x - 400, y' = y + 400
                vertices_coord3(0) = vertices_coord3(0) - 4000000;
                vertices_coord3(1) = vertices_coord3(1) + 4000000;
                vertices_coord3(2) = vertices_coord3(2) - 4000000;
                vertices_coord3(3) = vertices_coord3(3) + 4000000;
                vertices_coord3(4) = vertices_coord3(4) - 4000000;
                vertices_coord3(5) = vertices_coord3(5) + 4000000;
                vertices_coord3(6) = vertices_coord3(6) - 4000000;
                vertices_coord3(7) = vertices_coord3(7) + 4000000;
            }
        }//End of wrv_split == 1
        
        // fp = outputfile for wrv_split == 1
        if (wrv_split == 1) {
            fp = outputFile;
        }
        
        /*
         vertices_test[0] = vertices_coord(0);
         vertices_test[1] = vertices_coord(1);
         vertices_test[2] = vertices_coord(2);
         vertices_test[3] = vertices_coord(3);
         vertices_test[4] = vertices_coord(4);
         vertices_test[5] = vertices_coord(5);
         vertices_test[6] = vertices_coord(6);
         vertices_test[7] = vertices_coord(7);
         */
        
        piece = 1;
        pass = 0;
        
        //while loop to decide if it is ready to export
        while (piece == 1) {
            for (int q = 0; q< 4; q++) {
                vertices_ycoord(q) = vertices_coord3(q * 2 + 1);
                //vertices_ycoord_test[q] = vertices_coord3(q*2+1);					
            }
            
            // for the index
            std::ptrdiff_t i1, i2, i3, i4;
            
            //Sort out the max and min and avg
            num_1st = vertices_ycoord.maxCoeff(&i1);
            num_4th = vertices_ycoord.minCoeff(&i4);
            
            index_1st = i1;
            index_4th = i4;
            
            
            
            //Solve the situation first
            //3 pieces: 4 y coordinates are different
            //2 pieces: max=max or 
            //			min= min, then only 1 triangle and 1 rectangle
            //1 pieces: max=max and mint=min only 1 rectangle
            
            //Use mean to sort out 2nd and 3rd
            for (int q = 0; q< 4; q++) {
                vertices_ycoord(q) = vertices_ycoord(q) - num_4th;
                //vertices_ycoord_test[q] = vertices_ycoord(q);
            }
            
            //Read the number again for later comparison
            num_1st = vertices_ycoord.maxCoeff();
            num_4th = vertices_ycoord.minCoeff();
            
            //if there is 2 "0", then switch to the other way to sort out the order
            int test = 0;
            for (int q = 0; q< 4; q++) {
                if (num_4th == vertices_ycoord(q)) {
                    test = test + 1;
                }
            }
            
            if (test == 2) {
                // 2 min value = 0
                vertices_ycoord(index_1st) = 2 * num_1st;
                vertices_ycoord(index_4th) = 2 * num_1st;
                
                num_3rd = vertices_ycoord.minCoeff(&i3);
                index_3rd = i3;
                
                vertices_ycoord(index_3rd) = 2 * num_1st;
                
                num_2nd = vertices_ycoord.minCoeff(&i2);
                index_2nd = i2;
            }
            else {
                vertices_ycoord(index_1st) = 0;
                vertices_ycoord(index_4th) = 0;
                
                // The max in deltaY will be the 2nd, then 3rd.
                
                num_2nd = vertices_ycoord.maxCoeff(&i2);
                index_2nd = i2;
                vertices_ycoord(index_2nd) = 0;
                
                num_3rd = vertices_ycoord.maxCoeff(&i3);
                index_3rd = i3;
            }
            
            
            
            
            //Rearrange the vertices according to the y-axis, large to small
            vertices_coord2[0] = vertices_coord3((index_1st * 2 + 1) - 1);
            vertices_coord2[1] = vertices_coord3((index_1st * 2 + 1));
            vertices_coord2[2] = vertices_coord3((index_2nd * 2 + 1) - 1);
            vertices_coord2[3] = vertices_coord3(index_2nd * 2 + 1);
            vertices_coord2[4] = vertices_coord3((index_3rd * 2 + 1) - 1);
            vertices_coord2[5] = vertices_coord3(index_3rd * 2 + 1);
            vertices_coord2[6] = vertices_coord3((index_4th * 2 + 1) - 1);
            vertices_coord2[7] = vertices_coord3(index_4th * 2 + 1);
            
            //Check if it is satified the rule: trapezoid with 2 horizontal line or a triangle with horizontal line
            if (num_1st == num_2nd && num_3rd == num_4th) {
                //Trapezoid with 2 horizontal line or triangle with 1 horizontal line (Point upward and downward)
                pass = 1;
            }
            
            //Sort out if it is a trapezoid, triangle (point upward or downward)
            if (pass == 1) {
                
                //one_set = Sort(vertices_coord2, vertices_dose);
                
                //Smart_Cut2[6] = vertices_dose;
                
                if (vertices_coord2[0] < vertices_coord2[2]) {
                    if (vertices_coord2[4] < vertices_coord2[6]) {
                        //Trapezoid (LL LR UL UR = 3 4 1 2)
                        Smart_Cut2[0] = vertices_coord2[4];
                        Smart_Cut2[1] = vertices_coord2[5];
                        Smart_Cut2[2] = vertices_coord2[2];
                        Smart_Cut2[3] = vertices_coord2[3];
                        Smart_Cut2[4] = vertices_coord2[6];
                        Smart_Cut2[5] = vertices_coord2[0];
                        ////test_situation = 1;
                    }
                    else if (vertices_coord2[4] > vertices_coord2[6]) {
                        //Trapezoid (LL LR UL UR = 4 3 1 2)
                        Smart_Cut2[0] = vertices_coord2[6];
                        Smart_Cut2[1] = vertices_coord2[7];
                        Smart_Cut2[2] = vertices_coord2[2];
                        Smart_Cut2[3] = vertices_coord2[3];
                        Smart_Cut2[4] = vertices_coord2[4];
                        Smart_Cut2[5] = vertices_coord2[0];
                        //test_situation = 2;
                    }
                    else {
                        //vertices_coord2[4] = vertices_coord2[6]
                        //Triangle (LL UL UR = 4 1 2)
                        //Point downward
                        Smart_Cut2[0] = vertices_coord2[6];
                        Smart_Cut2[1] = vertices_coord2[7];
                        Smart_Cut2[2] = vertices_coord2[2];
                        Smart_Cut2[3] = vertices_coord2[3];
                        Smart_Cut2[4] = vertices_coord2[6];
                        Smart_Cut2[5] = vertices_coord2[0];
                        //test_situation = 3;
                    }
                }//if (vertices_coord2[0] < vertices_coord2[2])
                
                else if (vertices_coord2[0] > vertices_coord2[2]) {
                    if (vertices_coord2[4] < vertices_coord2[6]) {
                        //Trapezoid (LL LR UL UR = 3 4 2 1)
                        Smart_Cut2[0] = vertices_coord2[4];
                        Smart_Cut2[1] = vertices_coord2[5];
                        Smart_Cut2[2] = vertices_coord2[0];
                        Smart_Cut2[3] = vertices_coord2[1];
                        Smart_Cut2[4] = vertices_coord2[6];
                        Smart_Cut2[5] = vertices_coord2[2];
                        //test_situation = 4;
                    }
                    else if (vertices_coord2[4] > vertices_coord2[6]) {
                        //Trapezoid (LL LR UL UR = 4 3 2 1)
                        Smart_Cut2[0] = vertices_coord2[6];
                        Smart_Cut2[1] = vertices_coord2[7];
                        Smart_Cut2[2] = vertices_coord2[0];
                        Smart_Cut2[3] = vertices_coord2[1];
                        Smart_Cut2[4] = vertices_coord2[4];
                        Smart_Cut2[5] = vertices_coord2[2];
                        //test_situation = 5;
                    }
                    else {
                        //vertices_coord2[4] = vertices_coord2[6]
                        //Triangle (LL UL UR = 4 2 1)
                        //Point downward
                        Smart_Cut2[0] = vertices_coord2[6];
                        Smart_Cut2[1] = vertices_coord2[7];
                        Smart_Cut2[2] = vertices_coord2[0];
                        Smart_Cut2[3] = vertices_coord2[1];
                        Smart_Cut2[4] = vertices_coord2[6];
                        Smart_Cut2[5] = vertices_coord2[2];
                        //test_situation = 6;
                    }
                }//if (vertices_coord2[0] > vertices_coord2[2])
                
                else {//if (vertices_coord2[0] == vertices_coord2[2]){
                    //Triangle point upward
                    if (vertices_coord2[4] < vertices_coord2[6]) {
                        //Triangle (LL LR UR = 3 4 1)
                        Smart_Cut2[0] = vertices_coord2[4];
                        Smart_Cut2[1] = vertices_coord2[5];
                        Smart_Cut2[2] = vertices_coord2[0];
                        Smart_Cut2[3] = vertices_coord2[1];
                        Smart_Cut2[4] = vertices_coord2[6];
                        Smart_Cut2[5] = vertices_coord2[0];
                        //test_situation = 7;
                    }
                    else {
                        //Triangle (LL LR UR = 4 3 1)
                        Smart_Cut2[0] = vertices_coord2[6];
                        Smart_Cut2[1] = vertices_coord2[7];
                        Smart_Cut2[2] = vertices_coord2[0];
                        Smart_Cut2[3] = vertices_coord2[1];
                        Smart_Cut2[4] = vertices_coord2[4];
                        Smart_Cut2[5] = vertices_coord2[0];
                        //test_situation = 8;
                    }
                    //}//if (vertices_coord2[0] == vertices_coord2[2])
                }
                
                //SmartCut Check function
                //To make sure the output satify the definition.
                //y1 <= y3, x1 <= x2, x4 <= x3
                
                if (Smart_Cut2[0] > Smart_Cut2[4]) {
                    //x1 > x2
                    int32_t temp = Smart_Cut2[0];
                    Smart_Cut2[0] = Smart_Cut2[4];
                    Smart_Cut2[4] = temp;
                }
                
                if (Smart_Cut2[5] > Smart_Cut2[2]) {
                    //x4 > x3
                    int32_t temp = Smart_Cut2[2];
                    Smart_Cut2[2] = Smart_Cut2[5];
                    Smart_Cut2[5] = temp;
                }
                
                if (Smart_Cut2[1] > Smart_Cut2[3]) {
                    //y1 > x3
                    int32_t temp = Smart_Cut2[1];
                    Smart_Cut2[1] = Smart_Cut2[3];
                    Smart_Cut2[3] = temp;
                }
                
                
                //Save the one_set into SmartCut
                //SmartCut.conservativeResize((totalPoly+ 1)*7);
                
                
                int32_t test[6];
                
                for (int32_t m = 0; m < 6; m++) {
                    //double test = ((double) xArray[k*7+ m]/ (double)block_size)*block_size/5;
                    test[m] = (int32_t) ceil(((double)Smart_Cut2[m] / (double)block_size)*block_size / 10 / (double)pixel_size);
                }
                
                
                if (abs(test[1] - test[3]) > 1) {
                    
                    fprintf(fp, "trap/");
                    //dosage
                    fprintf(fp, "%ld ", vertices_dose);
                    
                    for (int m = 0; m < 6; m++) {
                        //double test = ((double) xArray[k*7+ m]/ (double)block_size)*block_size/5;
                        fprintf(fp, "%ld", (long) ceil(((double)Smart_Cut2[m] / (double)block_size)*block_size / 10 / (double)pixel_size));
                        if (m != 5)
                            fprintf(fp, " ");
                    }
                    //fprintf(fp, " %d", test_situation);
                    fprintf(fp, "\n");
                    //for (int q= 0; q< 7; q++){			
                    //	SmartCut(q+ totalPoly*7) = Smart_Cut2[q];
                    //	}
                    totalPoly = totalPoly + 1;
                }
                
                piece = 0;
                
            }//Pass =1
            else {
                // Pass ~= 1
                // Need to feed into the code to do the cut
                // Input: 8 vertices, output: 16 vertices
                // Always be on triangle and one trapezoid: 14 vertices+ 1 dose
                // Triangle always ready to be export! last 7 numbers, 6 vertices + 1 does
                
                //AfterCut = Cut(vertices_coord2, vertices_dose);
                
                
                for (int j = 0; j< 4; j++) {
                    vertices_coord2_y[j] = vertices_coord2[2 * j + 1];
                }
                
                //Already sort out by the order, thus if vertices_coord2_y[0] = vertices_coord2_y[1]
                //Then the 3rd point is the 2nd highest vertex and the triangle is at the bottom
                //If the vertices_coord2_y[0] ~= vertices_coord2_y[1], then the triangle is at the top
                // and the extra point share same y as 2nd point
                if (vertices_coord2_y[0] != vertices_coord2_y[1]) {
                    //Triangle is on the top, 2nd point y is used for the extra point
                    if (vertices_coord2[0] <= vertices_coord2[2]) {
                        // 2nd point is on the upper right
                        if (vertices_coord2[4] <= vertices_coord2[6]) {
                            if (vertices_coord2[0] != vertices_coord2[4]) {
                                //Use 1st and 3rd to find out the aL and bL for extra point
                                double aL = (double)(vertices_coord2[1] - vertices_coord2[5]) / (vertices_coord2[0] - vertices_coord2[4]);
                                double bL = (double)vertices_coord2[1] - aL*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[3] - bL) / aL;
                                extra_y = vertices_coord2[3];
                            }
                            else {
                                extra_x = vertices_coord2[4];
                                extra_y = vertices_coord2[3];
                            }
                            
                            Tri[0] = extra_x;
                            Tri[1] = extra_y;
                            Tri[2] = vertices_coord2[0];
                            Tri[3] = vertices_coord2[1];
                            Tri[4] = vertices_coord2[2];
                            Tri[5] = vertices_coord2[0];
                            //test_situation = 9;
                            
                            Trap[0] = vertices_coord2[4];
                            Trap[1] = vertices_coord2[5];
                            Trap[2] = vertices_coord2[6];
                            Trap[3] = vertices_coord2[7];
                            Trap[4] = vertices_coord2[2];
                            Trap[5] = vertices_coord2[3];
                            Trap[6] = extra_x;
                            Trap[7] = extra_y;
                        }
                        else {
                            //Use 1st and 4th to find out the aL and bL for extra point
                            if (vertices_coord2[0] != vertices_coord2[6]) {
                                
                                double aL = (double)(vertices_coord2[1] - vertices_coord2[7]) / (vertices_coord2[0] - vertices_coord2[6]);
                                double bL = (double)vertices_coord2[1] - aL*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[3] - bL) / aL;
                                extra_y = vertices_coord2[3];
                            }
                            else {
                                extra_x = vertices_coord2[6];
                                extra_y = vertices_coord2[3];
                            }
                            
                            Tri[0] = extra_x;
                            Tri[1] = extra_y;
                            Tri[2] = vertices_coord2[0];
                            Tri[3] = vertices_coord2[1];
                            Tri[4] = vertices_coord2[2];
                            Tri[5] = vertices_coord2[0];
                            //test_situation = 10;
                            
                            Trap[0] = vertices_coord2[6];
                            Trap[1] = vertices_coord2[7];
                            Trap[2] = vertices_coord2[4];
                            Trap[3] = vertices_coord2[5];
                            Trap[4] = vertices_coord2[2];
                            Trap[5] = vertices_coord2[3];
                            Trap[6] = extra_x;
                            Trap[7] = extra_y;
                        }
                    }//if(vertices_coord2[0] <= vertices_coord2[2])
                    else {
                        // 2nd point is on the upper left
                        if (vertices_coord2[4] <= vertices_coord2[6]) {
                            //Use 1st and 4th to find out the aR and bR for extra point
                            if (vertices_coord2[0] != vertices_coord2[6]) {
                                double aR = (double)(vertices_coord2[1] - vertices_coord2[7]) / (vertices_coord2[0] - vertices_coord2[6]);
                                double bR = (double)vertices_coord2[1] - aR*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[3] - bR) / aR;
                                extra_y = vertices_coord2[3];
                            }
                            else {
                                extra_x = vertices_coord2[0];
                                extra_y = vertices_coord2[3];
                            }
                            
                            
                            Tri[0] = vertices_coord2[2];
                            Tri[1] = vertices_coord2[3];
                            Tri[2] = vertices_coord2[0];
                            Tri[3] = vertices_coord2[1];
                            Tri[4] = extra_x;
                            Tri[5] = vertices_coord2[0];
                            //test_situation = 11;
                            
                            Trap[0] = vertices_coord2[4];
                            Trap[1] = vertices_coord2[5];
                            Trap[2] = vertices_coord2[6];
                            Trap[3] = vertices_coord2[7];
                            Trap[4] = extra_x;
                            Trap[5] = extra_y;
                            Trap[6] = vertices_coord2[2];
                            Trap[7] = vertices_coord2[3];
                            
                        }
                        else {
                            //Use 1st and 3rd to find out the aR and bR for extra point
                            if ((vertices_coord2[0] != vertices_coord2[4])) {
                                double aR = (double)(vertices_coord2[1] - vertices_coord2[5]) / (vertices_coord2[0] - vertices_coord2[4]);
                                double bR = (double)vertices_coord2[1] - aR*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[3] - bR) / aR;
                                extra_y = vertices_coord2[3];
                            }
                            else {
                                extra_x = vertices_coord2[0];
                                extra_y = vertices_coord2[3];
                            }
                            
                            Tri[0] = vertices_coord2[2];
                            Tri[1] = vertices_coord2[3];
                            Tri[2] = vertices_coord2[0];
                            Tri[3] = vertices_coord2[1];
                            Tri[4] = extra_x;
                            Tri[5] = vertices_coord2[0];
                            //test_situation = 12;
                            
                            Trap[0] = vertices_coord2[6];
                            Trap[1] = vertices_coord2[7];
                            Trap[2] = vertices_coord2[4];
                            Trap[3] = vertices_coord2[5];
                            Trap[4] = extra_x;
                            Trap[5] = extra_y;
                            Trap[6] = vertices_coord2[2];
                            Trap[7] = vertices_coord2[3];
                        }
                    }
                }//if (vertices_coord2_y[0] ~= vertices_coord2_y[1])-- Triangle is on the top, 2nd point y is used for the extra point
                else {
                    //Triangle is on the bottom, 3rd point y is used for the extra point
                    if (vertices_coord2[4] <= vertices_coord2[6]) {
                        //4th is the lower right
                        if (vertices_coord2[0] <= vertices_coord2[2]) {
                            // 2nd and the 4th point define aR and bR
                            if (vertices_coord2[2] != vertices_coord2[6]) {
                                double aR = (double)(vertices_coord2[3] - vertices_coord2[7]) / (vertices_coord2[2] - vertices_coord2[6]);
                                double bR = (double)vertices_coord2[3] - aR*vertices_coord2[2];
                                
                                extra_x = (vertices_coord2[5] - bR) / aR;
                                extra_y = vertices_coord2[5];
                            }
                            else {
                                extra_x = vertices_coord2[2];
                                extra_y = vertices_coord2[5];
                            }
                            
                            Tri[0] = vertices_coord2[6];
                            Tri[1] = vertices_coord2[7];
                            Tri[2] = extra_x;
                            Tri[3] = extra_y;
                            Tri[4] = vertices_coord2[6];
                            Tri[5] = vertices_coord2[4];
                            //test_situation = 13;
                            
                            Trap[0] = vertices_coord2[4];
                            Trap[1] = vertices_coord2[5];
                            Trap[2] = extra_x;
                            Trap[3] = extra_y;
                            Trap[4] = vertices_coord2[2];
                            Trap[5] = vertices_coord2[3];
                            Trap[6] = vertices_coord2[0];
                            Trap[7] = vertices_coord2[1];
                            
                        }
                        else {
                            //1st and 4th point define the aR and bR
                            if (vertices_coord2[0] != vertices_coord2[6]) {
                                double aR = (double)(vertices_coord2[1] - vertices_coord2[7]) / (vertices_coord2[0] - vertices_coord2[6]);
                                double bR = (double)vertices_coord2[1] - aR*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[5] - bR) / aR;
                                extra_y = vertices_coord2[5];
                            }
                            else {
                                extra_x = vertices_coord2[0];
                                extra_y = vertices_coord2[5];
                            }
                            
                            Tri[0] = vertices_coord2[6];
                            Tri[1] = vertices_coord2[7];
                            Tri[2] = extra_x;
                            Tri[3] = extra_y;
                            Tri[4] = vertices_coord2[6];
                            Tri[5] = vertices_coord2[4];
                            //test_situation = 14;
                            
                            Trap[0] = vertices_coord2[4];
                            Trap[1] = vertices_coord2[5];
                            Trap[2] = extra_x;
                            Trap[3] = extra_y;
                            Trap[4] = vertices_coord2[0];
                            Trap[5] = vertices_coord2[1];
                            Trap[6] = vertices_coord2[2];
                            Trap[7] = vertices_coord2[3];
                        }
                    }
                    else {
                        //4th is the lower left
                        if (vertices_coord2[0] <= vertices_coord2[2]) {
                            //1st and 4th define aL and bL
                            if (vertices_coord2[0] != vertices_coord2[6]) {
                                double aL = (double)(vertices_coord2[1] - vertices_coord2[7]) / (vertices_coord2[0] - vertices_coord2[6]);
                                double bL = (double)vertices_coord2[1] - aL*vertices_coord2[0];
                                
                                extra_x = (vertices_coord2[5] - bL) / aL;
                                extra_y = vertices_coord2[5];
                            }
                            else {
                                extra_x = vertices_coord2[0];
                                extra_y = vertices_coord2[5];
                            }
                            
                            Tri[0] = vertices_coord2[6];
                            Tri[1] = vertices_coord2[7];
                            Tri[2] = vertices_coord2[4];
                            Tri[3] = vertices_coord2[5];
                            Tri[4] = vertices_coord2[6];
                            Tri[5] = extra_x;
                            //test_situation = 15;
                            
                            Trap[0] = extra_x;
                            Trap[1] = extra_y;
                            Trap[2] = vertices_coord2[4];
                            Trap[3] = vertices_coord2[5];
                            Trap[4] = vertices_coord2[2];
                            Trap[5] = vertices_coord2[3];
                            Trap[6] = vertices_coord2[0];
                            Trap[7] = vertices_coord2[1];
                        }
                        else {
                            //2nd and 4th define the aL and bL
                            if (vertices_coord2[2] != vertices_coord2[6]) {
                                double aL = (double)(vertices_coord2[3] - vertices_coord2[7]) / (vertices_coord2[2] - vertices_coord2[6]);
                                double bL = (double)vertices_coord2[3] - aL*vertices_coord2[2];
                                
                                extra_x = (vertices_coord2[5] - bL) / aL;
                                extra_y = vertices_coord2[5];
                            }
                            else {
                                extra_x = vertices_coord2[2];
                                extra_y = vertices_coord2[5];
                            }
                            
                            Tri[0] = vertices_coord2[6];
                            Tri[1] = vertices_coord2[7];
                            Tri[2] = vertices_coord2[4];
                            Tri[3] = vertices_coord2[5];
                            Tri[4] = vertices_coord2[6];
                            Tri[5] = extra_x;
                            //test_situation = 16;
                            
                            Trap[0] = extra_x;
                            Trap[1] = extra_y;
                            Trap[2] = vertices_coord2[4];
                            Trap[3] = vertices_coord2[5];
                            Trap[4] = vertices_coord2[0];
                            Trap[5] = vertices_coord2[1];
                            Trap[6] = vertices_coord2[2];
                            Trap[7] = vertices_coord2[3];
                            
                        }
                    }
                }
                
                //one_set = Sort(vertices_coord2, vertices_dose);
                
                for (int q = 0; q< 8; q++) {
                    AfterCut[q] = Trap[q];
                }
                
                for (int q = 0; q< 6; q++) {
                    AfterCut[q + 8] = Tri[q];
                }
                
                //AfterCut[14] = vertices_dose;
                
                for (int q = 0; q< 8; q++) {
                    //Part has to be back to check if it is pass or not
                    vertices_coord3(q) = AfterCut[q];
                }
                
                //SmartCut Check function
                //To make sure the output satify the definition.
                //y1 <= y3, x1 <= x2, x4 <= x3
                
                if (AfterCut[8] > AfterCut[12]) {
                    //x1 > x2
                    int32_t temp = AfterCut[8];
                    AfterCut[8] = AfterCut[12];
                    AfterCut[12] = temp;
                }
                
                if (AfterCut[13] > AfterCut[10]) {
                    //x4 > x3
                    int32_t temp = AfterCut[10];
                    AfterCut[10] = AfterCut[13];
                    AfterCut[13] = temp;
                }
                
                if (AfterCut[9] >= AfterCut[11]) {
                    //y1 >= y3
                    int32_t temp = AfterCut[9];
                    AfterCut[9] = AfterCut[11];
                    AfterCut[11] = temp;
                }
                
                
                //Save the one_set into SmartCut
                
                //SmartCut.conservativeResize((totalPoly+ 1)*7);
                //for (int q= 0; q< 7; q++){			
                //	SmartCut(q+ totalPoly*7) = AfterCut[q+8];
                //}
                totalPoly = totalPoly + 1;
                
                int32_t test[6];
                for (int m = 0; m < 6; m++) {
                    //double test = ((double) xArray[k*7+ m]/ (double)block_size)*block_size/5;
                    test[m] = (int32_t) ceil(((double)AfterCut[m + 8] / (double)block_size)*block_size / 10);
                }
                
                if (abs(test[1] - test[3])> 1) {
                    
                    fprintf(fp, "trap/");
                    //dosage
                    fprintf(fp, "%ld ", vertices_dose);
                    
                    for (int m = 0; m < 6; m++) {
                        //double test = ((double) xArray[k*7+ m]/ (double)block_size)*block_size/5;
                        fprintf(fp, "%ld", (long)ceil(((double)AfterCut[m + 8] / (double)block_size)*block_size / 10 / (double)pixel_size));
                        if (m != 5)
                            fprintf(fp, " ");
                    }
                    //fprintf(fp, " %d", test_situation);
                    fprintf(fp, "\n");
                }
                
            }//End of else
        }//End of while loop
        
    }//for (long k = 0; k < numTrap; k++)
    
    if (first == 1) {
        //fclose(fp);
    }
    
    //Clear xArray memory
    
    //xArray = new double [1];
    
    
    //long long SmartCutSize = SmartCut.size();
    //long long* SmartCut_Info = new long long [SmartCutSize];
    
    //for (int index = 0; index < SmartCutSize; index++){
    //		SmartCut_Info[index] = SmartCut(index);
    //}
    
    //return SmartCut_Info;  	
}

// writes to text file for CAD_smartcut:
// Export as txt file for now 
// Has to switch between triangle and rectangular,
// which means from 6 vertices to 8 vertices
// To make sure the zoneplate is inside the block for wrv file
// The center of the zoneplate will shift to (block_size/2, block_size/2)

//Dummy function used by the older version
/*
 void writeToGDT_CAD(long* xArray, long totalpoly, double arrayC, char** argv_test, long zpID, double maxDose, double minDose, int block_size, int curl_on) {
 
	//    string dateStr  = "dateString";
	string fName = *argv_test;
	//    string cName    = "cellName";
	char    stringBuffer[100];
	FILE * fp;
	//long long numofTri;
	//long long numofRec;
 
 
	if ((fp = fopen(fName.c_str(), "wt")) == NULL)
 printf("Error");
	// Each trapezoid has 3 polygons
	// Index for dose
	int p = 0;
 
	fprintf(fp, "patdef 1000 %d %d %f %f 0 0\n", block_size / 10, block_size / 10, maxDose, minDose);
	fprintf(fp, "vepdef 20 %d %d\n", block_size / 10, block_size / 10);
	for (long k = 0; k < totalpoly; k++) {
 fprintf(fp, "trap/");
 //dosage
 fprintf(fp, "%ld ", xArray[k * 7 + 6]);
 
 for (long m = 0; m < 6; m++) {
 //double test = ((double) xArray[k*7+ m]/ (double)block_size)*block_size/5;
 fprintf(fp, "%ld", (long) floor(((double)xArray[k * 7 + m] / (double)block_size)*block_size / 10));
 if (m != arrayC - 1)
 fprintf(fp, " ");
 }
 fprintf(fp, "\n");
 
	}//for (long k = 0; k < totalpoly; k++)
 
	fclose(fp);
	//		printf("Total Rectangular in zone plate: %lld\n", numofRec);
	//		printf("Total Triangular in zone plate: %lld\n", numofTri);
 
 }
 */

void TrapTransform(float* ReproduceArc_array, int32_t &numTraps, int totalArc, double zTol, double bias, double & max_dose, double & min_dose, int File_format, int wrv_split,
                   double TrapDoesPower, char** argv_test, int32_t zpID, int block_size, int32_t & totalpoly, int curl_on) {
    // Seperate into 2 parts
    // Calculate the total number of Trap fist and build the array.
    // Then input the information
    // The seperation into different file tyep will happen here and write each arcs into file directly
    float r, theta, dr, dtheta, xc, yc, rs1, re1;
    //int name_size = NULL;
    //int* gds_Name_preamble = NULL;
    
    //For GDSII file format
    //unsigned char buffer2[4];
    unsigned char polyFormNBuffer[2]; // stores number of bytes
    //unsigned long numVertices;
    long numBytes;
    unsigned char gdsPost[8];
    unsigned char polyPre[16];
    unsigned char polyPost[4];
    unsigned char polyForm[2];
    //Only save each trpas
    //Only save each trpas & clock information
    
    float * trapCoords = new float[8];
    float * trapCoords_old = new float[8];
    int block_index;
    
    if (File_format == 3) {
        //WRV
        float * trapCoords = new float[9];
        float * trapCoords_old = new float[9];
        //Default is the center block
        block_index = 5;
    }
    
    int32_t totalSlices = 0;
    int32_t * trapsPerArc = new int32_t[totalArc];
    int32_t * cumulativeTrapsPerArc = new int32_t[totalArc];
    
    FILE* outputFile = NULL;
    FILE* outputFile_wrvsplit[9] = {};
    
    if (File_format == 1) {
        string fName = *argv_test;
        if ((outputFile = fopen(fName.c_str(), "wt")) == NULL)
            printf("Error");
    }
    else if (File_format == 2) {
        //GDSII
        
        char *gdsName;
        char *buffer;
        unsigned char sizeInfo[4];
        int Null_on = 0;
        
        int name_size = strlen(*argv_test) + 29;
        
        if (name_size % 2 != 0) {
            name_size = name_size + 1;
            Null_on = 1;
        }
        
        VectorXi gds_Name(name_size + 4);
        
        int gds_Name_length = gds_Name.size();
        
        gds_Name(0) = 0;
        gds_Name(1) = name_size + 4;
        gds_Name(2) = 6;
        gds_Name(3) = 6;
        
        if (Null_on == 0) {
            for (int q = 4; q < gds_Name_length; q++) {
                
                gds_Name(q) = *(*argv_test + 29 + (q - 4));
            }
        }
        else {
            for (int q = 4; q < gds_Name_length - 1; q++) {
                gds_Name(q) = *(*argv_test + 29 + (q - 4));
            }
            gds_Name(gds_Name_length - 1) = 0;
        }
        
        // Switch it back to string
        int* gds_Name_preamble = new int[name_size + 4];
        
        for (int j = 0; j < name_size + 4; j++) {
            *(gds_Name_preamble + j) = gds_Name(j);
        }
        
        
        strcat(*argv_test, ".gds");
        gdsName = *argv_test;
        
        if ((outputFile = fopen(gdsName, "wb")) == NULL)
            printf("Cannot open file.\n");
        
        // Define ambles:
        int gdspreamble[92] = { 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
            230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
            109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
            243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
            0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56 };
        //0, 11, 6, 6, 99, 110, 111, 110, 97, 109, 101};
        int gdspostamble[8] = { 0, 4, 7, 0, 0, 4, 4, 0 };
        int polypreamble[16] = { 0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0 };
        int polypostamble[4] = { 0, 4, 17, 0 };
        
        int polyBlockFormat[2] = { 16, 3 };
        
        // cast to char:
        unsigned char* gdsPre = new unsigned char[92 + name_size + 4];
        
        for (int k = 0; k < 92; k++)
            *(gdsPre + k) = (unsigned char)gdspreamble[k];
        for (int k = 92; k < 92 + name_size + 4; k++)
            *(gdsPre + k) = (unsigned char)gds_Name_preamble[k - 92];
        for (int k = 0; k < 8; k++)
            gdsPost[k] = (unsigned char)gdspostamble[k];
        for (int k = 0; k < 16; k++)
            polyPre[k] = (unsigned char)polypreamble[k];
        for (int k = 0; k < 4; k++)
            polyPost[k] = (unsigned char)polypostamble[k];
        for (int k = 0; k < 2; k++)
            polyForm[k] = (unsigned char)polyBlockFormat[k];
        
        
        // Write GDS preamble:
        for (int k = 0; k< 92 + name_size + 4; k++) {
            unsigned char word[1];
            word[0] = *(gdsPre + k);
            fwrite(word, sizeof(char), 1, outputFile);
        }
    }
    else if (File_format == 3) {
        //Determine if there is a need for extra blocks file
        if (wrv_split == 0) {
            strcat(*argv_test, ".wrv");
            string fName = *argv_test;
            if ((outputFile = fopen(fName.c_str(), "wt+")) == NULL)
                printf("Error");
        }
        else {
            // wrv_split == 1
            //Adding 9 wrv files with filenmae = filename"i".wrv
            //Top left block will be 1, bottom right block will be 9
            //The original/center block will be 5
            
            std::string fName_temp1(*argv_test);
            
            for (int k = 0; k < 9; k++) {
                //Put the numbers into the last element of fName_temp2
                std::string fName = fName_temp1;
                
                unsigned sz = fName.size();
                
                fName.resize(sz + 1, char('0' + k + 1));
                
                
                strcat(const_cast<char*>(fName.c_str()), ".wrv");
                if ((outputFile_wrvsplit[k] = fopen(fName.c_str(), "wt+")) == NULL) {
                    printf("Error");
                }
            }
        }
    }
    
    // Compute number of theta slices required in each zone		
    for (int index = 0; index < totalArc; index++) {
        
        r = ReproduceArc_array[0 + index * 6];
        theta = ReproduceArc_array[1 + index * 6];
        dr = ReproduceArc_array[2 + index * 6];
        dtheta = ReproduceArc_array[3 + index * 6];
        xc = ReproduceArc_array[4 + index * 6];
        yc = ReproduceArc_array[5 + index * 6];
        
        re1 = r + dr;
        
        // Keep track of number of trapezoids in each zone and cumlative:
        trapsPerArc[index] = (int32_t) ceil((dtheta / (2 * acos(-1 * (zTol*(dr) / re1 - 1)))));
        
        cumulativeTrapsPerArc[index] = totalSlices;
        totalSlices += trapsPerArc[index];
        
    }
    
    numTraps = totalSlices;
    //printf("Total traps in zone plate: %lld\n", totalSlices);
    //printf("Total traps in zone plate: %I64d\n", totalSlices);
    //Now we only save each traps and then pass it to next function
    //float * trapCoords = new float [totalSlices*8];
    
    VectorXf trapDose(totalSlices);
    float trapDose2;
    float trapClock;
    
    double s1x, s1y, e1x, e1y, s2x, s2y, e2x, e2y;
    int32_t numSlices;
    
    
    // Include a session for WRV only to find the max and min dose:
    
    if (File_format == 3) {
        
        for (int q = 0; q < totalArc; q++) {
            
            dr = ReproduceArc_array[2 + q * 6];
            numSlices = trapsPerArc[q];
            
            // Loop through each theta slice and compute coordinates
            for (int32_t k = 0; k < numSlices; k++) {
                
                // dr here is actually dr', which is dr' = dr- bias
                // dose equals to  dr/(dr - bias) = (dr'+ bias)/(dr'+ bias -bias)
                trapDose(cumulativeTrapsPerArc[q] + k) = (dr + bias / 1000) / (dr + bias / 1000 - bias / 1000);
            }//End of k
            
        }//End of q
        
        //This part is used to include trapDose_power into the consideration
        // we have to switch between array and matrix in order to do the power calculation.
        trapDose = ((trapDose.array()).pow(TrapDoesPower)).matrix();
        
        //Normalized the dose and scale it up to 2^16
        max_dose = trapDose.maxCoeff();
        min_dose = trapDose.minCoeff();
    }//end of if (File_format == 3){
    
    double first = 0;
    double first_wrvSplit[9] = {};
    
    // Loop through each zone
    for (int q = 0; q < totalArc; q++) {
        
        
        
        numSlices = trapsPerArc[q];
        r = ReproduceArc_array[0 + q * 6];
        theta = ReproduceArc_array[1 + q * 6];
        dr = ReproduceArc_array[2 + q * 6];
        dtheta = ReproduceArc_array[3 + q * 6];
        xc = ReproduceArc_array[4 + q * 6];
        yc = ReproduceArc_array[5 + q * 6];
        re1 = r + dr;
        
        if (File_format == 3) {
            
            trapDose2 = (dr + bias / 1000) / (dr + bias / 1000 - bias / 1000);
            float a = 1 / trapDose2;
            double b = 1 / max_dose;
            double c = 1 / min_dose;
            // In case we don't have bias....
            if (c != b) {
                trapClock = ((a - b) / (c - b))*(pow((long double)2, 16) - 1);
                //To deal with the situation that floor(trapClock) = -1;
                if (trapClock < 0) {
                    trapClock = 0;
                }
            }
            else {
                trapClock = pow((long double)2, 16) - 1;
            }
            
            trapCoords_old[8] = trapClock;
            trapCoords[8] = trapClock;
        }
        
        // Loop through each theta slice and compute coordinates
        for (int32_t k = 0; k < numSlices; k++) {
            // theta- dtheta/2: to start from the edge, check arc definition
            
            double trap_theta = theta - dtheta / 2 + k* dtheta / ((double)numSlices);
            
            // Compute zone locations at this theta
            
            s1x = r*cos(trap_theta) + xc;
            s1y = r*sin(trap_theta) + yc;
            e1x = re1*cos(trap_theta) + xc;
            e1y = re1*sin(trap_theta) + yc;
            
            trapCoords[0] = s1x;
            trapCoords[1] = s1y;
            trapCoords[6] = e1x;
            trapCoords[7] = e1y;
            
            if (k > 0) {
                trapCoords_old[4] = trapCoords[6];
                trapCoords_old[5] = trapCoords[7];
                trapCoords_old[2] = trapCoords[0];
                trapCoords_old[3] = trapCoords[1];
            }
            if (k == (numSlices - 1)) {
                
                trapCoords[4] = re1*cos(theta + dtheta / 2) + xc;
                trapCoords[5] = re1*sin(theta + dtheta / 2) + yc;
                trapCoords[2] = r*cos(theta + dtheta / 2) + xc;
                trapCoords[3] = r*sin(theta + dtheta / 2) + yc;
            }
            // Pass it to the following function here once you have all the trapCoords for each one:
            
            if (q == (totalArc - 1) && k == (numSlices - 1)) {
                //Vary last piece
                first = 1;
            }
            //if (k > 0 && k != (numSlices - 1)) {
            if (k > 0 ) {
                if (File_format == 1) {
                    // txt file with only coordiantes
                    // writes trapezoid coordinates to a text file:
                    writeToGDSCell_Trap(trapCoords_old, 1, 8, argv_test, zpID, first, outputFile, curl_on);
                }
                else if (File_format == 2) {
                    // GDSII file
                    //Switch to a function to write in the form of gds directly
                    //mexFunction(trapCoords_old, 1, 8, argv_test, zpID, first, outputFile, name_size);
                    for (int k = 0; k < 1; k++) {
                        
                        // write polypreamble
                        fwrite(polyPre, sizeof(char), 16, outputFile);
                        
                        // write polyform:
                        encode16(44, polyFormNBuffer);
                        //double a = decode32(polyFormNBuffer);
                        
                        //printf("numBytes: [%d]\n", (int)numBytes);
                        //printf("polyFormNBuffer: [%d, %d]\n", (int)polyFormNBuffer[0], (int)polyFormNBuffer[1]);
                        
                        fwrite(polyFormNBuffer, sizeof(char), 2, outputFile);
                        fwrite(polyForm, sizeof(char), 2, outputFile);
                        
                        // read and write block
                        unsigned char dataBuffer2[40];
                        
                        for (int m = 0; m < 8; m++) {
                            int32_t num = trapCoords_old[1 * m + k] * 10000;
                            unsigned char dataBuffer[4];
                            encode32(num, dataBuffer);
                            for (int q = 0; q < 4; q++) {
                                dataBuffer2[m * 4 + q] = dataBuffer[q];
                            }
                        }
                        
                        //Rewrite the 1st verticies
                        for (int m = 0; m < 2; m++) {
                            int32_t num = trapCoords_old[1 * m + k] * 10000;
                            unsigned char dataBuffer[4];
                            encode32(num, dataBuffer);
                            for (int q = 0; q< 4; q++) {
                                dataBuffer2[(m + 8) * 4 + q] = dataBuffer[q];
                            }
                        }
                        
                        fwrite(dataBuffer2, 1, 40, outputFile);
                        
                        //thisPolyBlock += polypostamble;
                        fwrite(polyPost, sizeof(char), 4, outputFile);
                    }
                }
                else if (File_format == 3) {
                    // CAD file
                    //Do smart cut for CAD tool
                    //CAD_smartcut = CAD_tool(trapCoords, numTraps, 8, argv_test, zpID, block_size, totalpoly, argv_test, zpID, maxDoes, minDoes);
                    
                    block_index = BlockNumberFinder(trapCoords_old, block_size);
                    CAD_tool(trapCoords_old, 1, 8, argv_test, zpID, block_size, wrv_split, block_index, totalpoly, max_dose, min_dose, first, first_wrvSplit
                             , outputFile, outputFile_wrvsplit);
                }
                
                //initially: first = 0
                //After: first = 2
                //Last: first = 1
                if (wrv_split == 0) {
                    if (first != 1) {
                        first = 2;
                    }
                }
                else {
                    //if wrv split
                    if (first != 1) {
                        first = 2;
                    }
                    
                    if (first_wrvSplit[block_index - 1] != 1) {
                        first_wrvSplit[block_index - 1] = 2;
                    }
                }
            }//End of file input k > 0
            
            if (k == (numSlices - 1)) {
                if (File_format == 1) {
                    // txt file with only coordiantes
                    // writes trapezoid coordinates to a text file:
                    writeToGDSCell_Trap(trapCoords, 1, 8, argv_test, zpID, first, outputFile, curl_on);
                }
                else if (File_format == 2) {
                    // GDSII file
                    //Switch to a function to write in the form of gds directly
                    //mexFunction(trapCoords, 1, 8, argv_test, zpID, first, outputFile, name_size);
                    for (int k = 0; k < 1; k++) {
                        
                        // write polypreamble
                        fwrite(polyPre, sizeof(char), 16, outputFile);
                        
                        // write polyform:
                        encode16(44, polyFormNBuffer);
                        //double a = decode32(polyFormNBuffer);
                        
                        //printf("numBytes: [%d]\n", (int)numBytes);
                        //printf("polyFormNBuffer: [%d, %d]\n", (int)polyFormNBuffer[0], (int)polyFormNBuffer[1]);
                        
                        fwrite(polyFormNBuffer, sizeof(char), 2, outputFile);
                        fwrite(polyForm, sizeof(char), 2, outputFile);
                        
                        // read and write block
                        unsigned char dataBuffer2[40];
                        
                        for (int m = 0; m < 8; m++) {
                            int32_t num = trapCoords[1 * m + k] * 10000;
                            unsigned char dataBuffer[4];
                            encode32(num, dataBuffer);
                            for (int q = 0; q < 4; q++) {
                                dataBuffer2[m * 4 + q] = dataBuffer[q];
                            }
                        }
                        
                        //Rewrite the 1st verticies
                        for (int m = 0; m < 2; m++) {
                            int32_t num = trapCoords[1 * m + k] * 10000;
                            unsigned char dataBuffer[4];
                            encode32(num, dataBuffer);
                            for (int q = 0; q< 4; q++) {
                                dataBuffer2[(m + 8) * 4 + q] = dataBuffer[q];
                            }
                        }
                        
                        fwrite(dataBuffer2, 1, 40, outputFile);
                        
                        //thisPolyBlock += polypostamble;
                        fwrite(polyPost, sizeof(char), 4, outputFile);
                    }
                }
                else if (File_format == 3) {
                    // CAD file
                    //Do smart cut for CAD tool
                    //CAD_smartcut = CAD_tool(trapCoords, numTraps, 8, argv_test, zpID, block_size, totalpoly, argv_test, zpID, maxDoes, minDoes);
                    block_index = BlockNumberFinder(trapCoords, block_size);
                    CAD_tool(trapCoords, 1, 8, argv_test, zpID, block_size, wrv_split, block_index, totalpoly, max_dose, min_dose, first, first_wrvSplit,
                             outputFile, outputFile_wrvsplit);
                }
                
                //initially: first = 0
                //After: first = 2
                //Last: first = 1
                if (wrv_split == 0) {
                    if (first != 1) {
                        first = 2;
                    }
                }
                else {
                    //If wrv_split == 1
                    if (first != 1) {
                        first = 2;
                    }
                    if (first_wrvSplit[block_index - 1] != 1) {
                        first_wrvSplit[block_index - 1] = 2;
                    }
                }
            }//End of file input (k == numSlices - 1)
            
            if (first != 1) {
                trapCoords_old[0] = s1x;
                trapCoords_old[1] = s1y;
                trapCoords_old[6] = e1x;
                trapCoords_old[7] = e1y;
            }
        }//End of k
        char    stringBuffer[100];
        
        if (curl_on == 1) {
            if (q % 10 == 0) {
                sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                        zpID, (((float)(q+1) / (float)totalArc) / 2 + 0.5));
                system(stringBuffer);
                //printf("Counting m0.5: %0.3f \n", (((float)(q + 1) / (float)totalArc) / 2 + 0.5));
                
            } 
        }
        
    }//End of q
    
    //delete[] trapCoords_old;
    //delete[] trapCoords;
    
    /*if (first == 1 && File_format == 2) {
     // Write GDS postamble:
     fwrite(gdsPost, sizeof(char), 8, outputFile);
     fclose(outputFile);		
     }
     if (wrv_split == 1) {
     //You have to close all 9 files
     for (int k = 0; k < 9; k++) {
     fclose(outputFile_wrvsplit[k]);
     }
     }
     else {
     fclose(outputFile);
     }
     */
    
    //The old version that is work without crash...
    
    if (first == 1 && File_format == 2) {
        // Write GDS postamble:
        fwrite(gdsPost, sizeof(char), 8, outputFile);
        //fclose(op);
    }
    if (wrv_split == 1) {
        //You have to close all 9 files
        for (int k = 0; k < 9; k++) {
            fclose(outputFile_wrvsplit[k]);
        }
        char stringBuffer[100];
        if (curl_on == 1) {
            sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                    zpID, (float) 1.0);
            system(stringBuffer);
            //printf("Done: 1 \n");
        }
    }
    else {
        fclose(outputFile);
        char stringBuffer[100];
        if (curl_on == 1) {
            sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
                    zpID, (float) 1.0);
            system(stringBuffer);
            //printf("Done: 1 \n");
        }
    }
    
    
    //Provide fileformat information
    if (File_format == 1) {
        printf("Save as trapezoid information in txt file\n");
    }
    else if (File_format == 2) {
        printf("Save as GDSII file \n");
    }
    else if (File_format == 3) {
        printf("Save as wrv file \n");
    }
}



int main(int argc, char** argv)
{
    /*
     
     p : object distance
     q : image distance
     f : focal length
     M: magnification
     
     p =  (k+1)f;
     q =  ((k+1)/k)*f;
     f = (p*q)/(p+q);
     k = q/p;
     double f_mm           = 4.44;
     p & q in um
     double mag            = 10000;
     */
    
    /*
     double zTol             = .001;
     double lambda_nm        = 1.77;
     double p                =  10001*4440;
     double q                = (10001/10000)*4440;
     double obscurationSigma = 0.01;
     double NA_P               = 0.005;
     double orders[2]        = {4, 1};
     double alpha            = 0;
     double beta             = M_PI/4;
     double offaxis_R1       = 0;
     double offaxis_R2       = 2;
     int offaxis_numpoints   = 1000;
     double offaxis_theta    = 2*M_PI;
     */
    srand((unsigned)time(NULL));
    
    if (argc == 1) {
        
        printf("No input parameter!!! \n");
        system("PAUSE");
        return 0;
        
    }
    
    char** argv_test;
    
    argv_test = argv + 1;
    
    double zTol = atof(*(argv_test++));
    double lambda_nm = atof(*(argv_test++));
    double p = atof(*(argv_test++));
    double q = atof(*(argv_test++));
    double obscurationSigma = atof(*(argv_test++));
    double NA_P = atof(*(argv_test++));
    double number_of_orders = atof(*(argv_test++));
    
    VectorXd orders_eigen(1);
    orders_eigen.conservativeResize(2 * number_of_orders);
    for (int k = 0; k < number_of_orders; k++) {
        
        orders_eigen(k) = atof(*(argv_test++));
        orders_eigen(k + number_of_orders) = atof(*(argv_test++));
    }
    
    int order_size = orders_eigen.size();
    double* order = new  double[order_size];
    
    for (int k = 0; k < order_size; k++) {
        
        order[k] = orders_eigen(k);
    }
    
    // Need another array to finish transport,except using vectorXd
    // gamma: the angle rotate on the plane perpendicular to optica axis
    
    double alpha = atof(*(argv_test++));
    double beta = atof(*(argv_test++));
    double gamma = atof(*(argv_test++));
    
    double CRA = atof(*(argv_test++));
    double NA_C = atof(*(argv_test++));
    // NoP: number of parts (total)
    // IoP: index of parts (interested)
    double PC_term = atof(*(argv_test++));
    double APD = atof(*(argv_test++));
    double APD_window = atof(*(argv_test++));
    double InnerSigma = atof(*(argv_test++));
    double OuterSigma = atof(*(argv_test++));
    double bias_nm = atof(*(argv_test++));
    int File_format = atof(*(argv_test++));
    int Opposite_Tone = atof(*(argv_test++));
    int Free_standing = atof(*(argv_test++));
    double W = atof(*(argv_test++));
    double T = atof(*(argv_test++));
    double TrapDoesPower = atof(*(argv_test++));
    //1 million pixel with a pixel size 500 pm (0.5 nm)
    // Thus it transfer to 500 um block size
    int block_size = atof(*(argv_test++));
    double NoP = atol(*(argv_test++));
    double IoP = atol(*(argv_test++));
    int32_t zpID = atol(*(argv_test++));
    int curl_on = atol(*(argv_test++));
    
    
    //double f;       // focus in microns    
    double lambda;  // lambda in microns    
    double NA_Pt;     // tan(theta)
    double D;       // diameter of zone plate
    double R0;      // radius of zone plate
    double r0;      // radius of obscuration
    double ri;      // Inner radius of PC/APD
    double ro;      // Outer radius of PC/APD
    
    
    int ns;
    int ne;
    int ni = 0;
    int no = 0;
    int offaxis_numpoints = 100000;
    int wrv_split = 0;
    double x_test, y_test, R_offZP;
    
    double Rc = p*tan((CRA / 180)*M_PI);
    double offaxis_R2 = p / 2 * ((sin((CRA / 180.0)*M_PI) + NA_C) / sqrt(1 - (sin((CRA / 180.0)*M_PI) + NA_C)*(sin((CRA / 180.0)*M_PI) + NA_C)) - (sin((CRA / 180.0)*M_PI) - NA_C) / sqrt(1 - (sin((CRA / 180.0)*M_PI) - NA_C)*(sin((CRA / 180.0)*M_PI) - NA_C)));
    //double Roff = p*(tan((CRA / 180.0)*M_PI)- tan((1 / 180.0)*M_PI));
    
    x_test = (p*(NA_C*sin(0) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(0) - NA_C*NA_C)));
    y_test = (p*((sin((CRA / 180)*M_PI) + NA_C*cos(0)) / sqrt((cos((CRA / 180)*M_PI)*cos((CRA / 180)*M_PI)) - 2 * NA_C*sin((CRA / 180)*M_PI)*cos(0) - NA_C*NA_C)));
    R_offZP = sqrt(x_test*x_test + (y_test - Rc)* (y_test - Rc));
    
    double Doff = 2 * R_offZP;
    
    // convert to microns:\
    
    //f = f_mm*1000;
    // lambda in um;
    lambda = lambda_nm / 1000;
    
    NA_Pt = tan(asin(NA_P));
    
    D = 2 * (p*q / (p + q))*tan(asin(NA_P));;
    
    R0 = D / 2;
    
    //From block_size x resolution limit (1000pm = 1nm)
    block_size = block_size * 10;
    
    // Define the ZP starting point, and also the PC/APD region
    // Have to differentiate case for on/off-axis
    r0 = R0* obscurationSigma;
    ri = R0* InnerSigma;
    ro = R0* OuterSigma;
    
    
    ns = nfinder(r0, alpha, beta, gamma, lambda, p, q);
    ne = nfinder(R0, alpha, beta, gamma, lambda, p, q);
    
    //01232017
    //If you turn off this part CRA == 0
    //You can extend the off-axis feature to the zoneplate
    //with Phase contrast of apodization in their central region
    //
    
    if (CRA == 0) {
        if (OuterSigma != 0) {
            no = nfinder(ro, alpha, beta, gamma, lambda, p, q);
        }
        
        if (InnerSigma != 0) {
            ni = nfinder(ri, alpha, beta, gamma, lambda, p, q);
        }
    }
    
    
    if (ns <= 0)
        ns = 1;
    
    //In case there isn't any ZP to build
    
    if (ns >= ne) {
        printf("\nNs > Ne, no ZPs has to be built!!!\n");
        system("PAUSE");
        return 0;
    }
    if (no != 0) {
        
        if (ni >= no) {
            printf("\nOuter Sigma has to be larger than Inner Sigma!!!\n");
            system("PAUSE");
            return 0;
        }
    }
    
    if (NoP == 1 && IoP != 1) {
        printf("\nIndex of part has to be 1 if number of parts is 1 !\n");
        system("PAUSE");
        return 0;
        
    }
    
    //For Free Standing mode 2 only: Only print the smaller piece
    // Switch W and T
    
    if (Free_standing == 2) {
        double W_temp = W;
        W = T;
        T = W_temp;
    }
    
    // For Free Standing case 3: print the opposite tone.
    if (Opposite_Tone == 1) {
        ns = ns - 1;
        ne = ne - 1;
        if (ns == 0) {
            ns = 2;
        }
        if (ni != 0 && no != 0) {
            ni = ni - 1;
            no = no - 1;
            if (ni == 0) {
                ni = 2;
            }
        }
        
    }
    
    printf("\nCreating zone plate:\nNA_P \t\t\t= %0.3f\n", NA_P);
    printf("NA_child \t\t= %f \n", NA_C);
    printf("lambda \t\t\t= %0.3f nm\n", lambda_nm);
    printf("objective length \t= %g um\n", p);
    printf("image length \t\t= %g um\n", q);
    
    //Aberration order
    for (int j = 0; j < order_size / 2; j++) {
        printf("Zernike Number \t\t= %f \n", order[j]);
        printf("Weight \t\t\t= %f waves\n", order[j + order_size / 2]);
    }
    
    printf("titlted angle alpha \t= %g in radian\n", alpha);
    printf("titlted angle beta \t= %g in radian\n", beta);
    printf("titlted angle gamma \t= %g in radian\n", gamma);
    printf("Max zone phase error\t= lambda/%d\n", (int)(1 / zTol));
    printf("Starting Zone N \t= %d\n", ns);
    printf("Ending Zone N \t\t= %d\n", ne);
    if (CRA == 0) {
        printf("PCZP Starting Zone N \t= %d\n", ni);
        printf("PCZP Ending Zone N \t= %d\n", no);
    }
    printf("Zone plate diameter \t= %0.2f um\n", D);
    
    if (CRA != 0) {
        printf("CRA \t\t\t= %f degree \n", CRA);
        printf("Rc\t\t\t= %g um \n", p*tan(CRA / 180 * M_PI));
        printf("Radius of off-axis pattern = %g um\n", offaxis_R2);
        printf("Diameter of off-axis pattern = %g um\n", Doff);
        printf("Radius of PC/APD region = %g um\n", offaxis_R2* OuterSigma);
    }
    
    printf("Extra Phase shift \t= %f in Radian\n", PC_term);
    printf("Apodized term \t\t= %f %% transmission\n", APD * 100);
    //If we initiate specific apodization window
    if (APD_window != 0) {
        printf("Apodized window \t= %f \n", APD_window);
    }
    printf("Bias term \t\t= %f nm\n", bias_nm);
    printf("Free Standing \t\t= %d \n", Free_standing);
    printf("W\t\t\t= %f dr \n", W);
    printf("T\t\t\t= %f dr \n", T);
    printf("Trap Dose Power \t= %f\n", TrapDoesPower);
    printf("Block sieze for wrv file: %d\n", block_size);
    printf("Generate %f part out of total %f parts.\n", IoP, NoP);
    //printf("Transfer to Trapozoid (Yes:1, No: 0) = %d\n", Trap_on);
    if (Opposite_Tone == 1) {
        printf("Print in opposie tone!");
    }
    // Over the limit of the block size for WRV file
    // Intiate the block split into 9 different block
    if (D >= 400 && File_format == 3 && CRA == 0) {
        wrv_split = 1;
        printf("On-axis Zoneplate is larger than WRV block size, split into 9 different blocks/files");
    }
    else if (Doff >= 400 && File_format == 3 && CRA != 0) {
        wrv_split = 1;
        printf("Off-axis Zoneplate is larger than WRV block size, split into 9 different blocks/files");
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
    
    printf("\n");
    
    int rowLen = number_of_orders;
    
    // For titling
    double* CT_array;
    double* CT_array2; // For PC/APD case
    double* CT_array3; // For PC/APD case (InnerSigma)
    // For off-axis situation
    double* nandtheta_array;
    double* nandtheta_array2;// For PC/APD case
    double* nandtheta_array3;// For PC/APD case (InnerSigma)
    // For the final ZP information
    float* ReproduceArc_array;
    // For the final ZP information in trapozoid
    float * trapCoords = NULL;
    //long long * CAD_smartcut;
    
    
    
    int totalArcs = 0;
    int32_t numTraps = 0;
    int32_t totalpoly = 0;
    int CT_array_size = 0;
    int nandtheta_array_size = 0;
    int CT_array_size2 = 0;
    int nandtheta_array_size2 = 0;
    int CT_array_size3 = 0;
    int nandtheta_array_size3 = 0;
    double offaxis_R_Outer = 0;
    double offaxis_R_Inner = 0;
    int offaxis_ZoneS = 0;
    int offaxis_ZoneE = 0;
    int offaxis_aberration = 0;
    double maxDoes;
    double minDoes;
    
    
    //int numZones = 0;
    //double*data_array;
    
    /* Test only
     double r_test = R0*0.5;
     double n_test = nfinder(r_test, alpha, beta, lambda, p, q);
     double x_test1  = sqrt(n_test * lambda * (p*q)/(p+q));
     double rs = secantSolve(x_test1, 0, orders, 1, n_test, p, lambda, q, alpha, beta, R0, 1, 20, 0.01, 1.57, 0);
     double x_test2  = sqrt((n_test+1) * lambda * (p*q)/(p+q));
     double re = secantSolve(x_test2, 0, orders, 1, n_test+1, p, lambda, q, alpha, beta, R0, 1, 20, 0.01, 1.57, 0);
     double r_final = re-rs;
     */
    
    // on-Axis case
    if (CRA == 0) {
        ReproduceArc_array = HWcustomZP(zTol, lambda, p, q, alpha, beta, gamma, ns, ni, no, ne, R0, order, rowLen, totalArcs, PC_term, APD, bias_nm, zpID, offaxis_aberration, CRA, NA_C, File_format, Free_standing, W, T, argv_test, NoP, IoP, curl_on);
    }
    //off-axis case with PC/APD/APD window
    else if (PC_term != 0 || APD != 1 || APD_window != 0) {
        
        //old off-axis function, temporarily suspend Henry Wang 01/20/2017
        // off-Axis case
        // No PC & aberration
        if (number_of_orders == 0) {
            CT_array = CoordinateTransform(0, 0, 0, 0, offaxis_numpoints, CT_array_size, PC_term, 1, NA_C, p, q, CRA, APD);
            nandtheta_array = define_function(CT_array, p, q, lambda, NA_P, 0, 0, 0, R0, order, rowLen, CT_array_size, nandtheta_array_size, offaxis_ZoneS, offaxis_ZoneE, offaxis_aberration, CRA, NA_C);
            nandtheta_array2 = new double[0];
            nandtheta_array3 = new double[0];
        }
        else {
            offaxis_aberration = 1;
            CT_array = CoordinateTransform(0, 0, 0, 0, offaxis_numpoints, CT_array_size, PC_term, 1, NA_C, p, q, CRA, APD);
            nandtheta_array = define_function(CT_array, p, q, lambda, NA_P, 0, 0, 0, R0, order, rowLen, CT_array_size, nandtheta_array_size, offaxis_ZoneS, offaxis_ZoneE, offaxis_aberration, CRA, NA_C);
            nandtheta_array2 = new double[0];
            nandtheta_array3 = new double[0];
            
        }
        
        // With PC or APD  
        if (PC_term != 0 || APD != 1) {
            // Disk illumination
            if (OuterSigma != 0) {
                
                offaxis_R_Outer = offaxis_R2* OuterSigma;
                CT_array2 = CoordinateTransform(offaxis_R_Outer, alpha, beta, gamma, offaxis_numpoints, CT_array_size2, PC_term, 0, NA_C, p, q, CRA, APD);
                nandtheta_array2 = define_function(CT_array2, p, q, lambda, NA_P, alpha, beta, gamma, R0, order, rowLen, CT_array_size2, nandtheta_array_size2, offaxis_ZoneS, offaxis_ZoneE, offaxis_aberration, CRA, NA_C);
                nandtheta_array3 = new double[0];
                
            }
            // Annular illumination
            if (InnerSigma != 0) {
                
                double offaxis_R_Inner = offaxis_R2* InnerSigma;
                CT_array3 = CoordinateTransform(offaxis_R_Inner, alpha, beta, gamma, offaxis_numpoints, CT_array_size3, PC_term, 0, NA_C, p, q, CRA, APD);
                nandtheta_array3 = define_function(CT_array3, p, q, lambda, NA_P, alpha, beta, gamma, R0, order, rowLen, CT_array_size3, nandtheta_array_size3, offaxis_ZoneS, offaxis_ZoneE, offaxis_aberration, CRA, NA_C);
                
            }
        }
        
        ReproduceArc_array = Control(nandtheta_array, nandtheta_array2, nandtheta_array3, lambda, p, q, alpha, beta, gamma, ns, ne, obscurationSigma,
                                     InnerSigma, OuterSigma, NA_P, order, rowLen, totalArcs, zTol, nandtheta_array_size, nandtheta_array_size2,
                                     nandtheta_array_size3, PC_term, APD, bias_nm, zpID, offaxis_aberration, CRA, NA_C, File_format, argv_test, NoP, IoP, APD_window, curl_on);
        
    }
    else{
        //off-axis case with aberration or free standing only
        
        if (number_of_orders != 0) {
            offaxis_aberration = 1;
        }
        ReproduceArc_array = HWcustomZP_offaxis(zTol, lambda, p, q, alpha, beta, gamma, ns, ni, no, ne, R0, order, rowLen, totalArcs, PC_term, APD, bias_nm, zpID, offaxis_aberration, CRA, NA_C, File_format, Free_standing, W, T, argv_test, NoP, IoP, curl_on);
    }
    
    
    
    printf("Zone plate computation complete\n");
    //Status bar should reach 50% at this point
    
    //For the file format, Arc txt GDSII WRV file will be name as File_format = 0, 1, 2, 3
    if (File_format == 0) {
        // Arc
        //   
        writeToGDSCell(ReproduceArc_array, totalArcs, 6, argv_test, zpID, curl_on);
        printf("Save as Arc information\n");
    }
    else if (File_format != 0 && File_format != 4) {
        
        //trapCoords = TrapTransform(ReproduceArc_array, numTraps, totalArcs, zTol, bias_nm, maxDoes, minDoes, File_format, argv_test, zpID);
        TrapTransform(ReproduceArc_array, numTraps, totalArcs, zTol, bias_nm, maxDoes, minDoes, File_format, wrv_split, TrapDoesPower, argv_test, zpID,
                      block_size, totalpoly, curl_on);
    }
    else {
        printf("Save as Zone Number and Clock information\n");
    }
    
    //char stringBuffer[100];
    //if (curl_on == 1) {
    //sprintf(stringBuffer, "curl \"http://joanie2.msd.lbl.gov/zpdev/index.php?r=zpStatus/putProgress&zpID=%ld&progress=%0.3f\"",
    //	zpID,(float) 1);
    //printf("Done: 1 \n");
    //system(stringBuffer);
    //}
    
    system("PAUSE");
    return 0;
    
}
