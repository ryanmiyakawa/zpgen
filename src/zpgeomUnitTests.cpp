#include <iostream>
#include <cmath>
#include "zpUtils.h"

int main() {
    const double lambda_nm = 13.5;
    const double na = 0.33;
    const double T_MIN_nm = lambda_nm / na;
    const double T_MIN_um = T_MIN_nm / 1000;
    const double lambda_um = lambda_nm / 1000;
    double beta, n[3], u[3], r[3], f[2], k[2], bx[3], by[3], bz[3], U[3], r2[3];
    double fx, fy, rx, p = 1e3, q = 1e4;
    double *fxy = new double[2];

    // ======== 1.1 UNIT TEST nominal zone plate
    beta = 0.0;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    u[0] = 0; u[1] = 0; u[2] = p;
    fx = 1 / T_MIN_um * 0.5;
    fy = 0;
    fxy[0] = fx; fxy[1] = fy;
    freq2zpXYZ(fxy, n, u, lambda_um, r);
    rx = tan(asin(lambda_nm / T_MIN_nm / 2));
    std::cout << "UNIT TEST 1.1 (nominal zp): Geometric: " << rx << ", fn1: " << r[0] / 1000 << std::endl;

    // UNIT TEST 1.2 (tilt-x zp)
    beta = 0.1;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    freq2zpXYZ(fxy, n, u, lambda_um, r);
    printf("UNIT TEST 1.2 (tilt-x zp): Geometric: %0.3f, fn1: %0.3f\n", rx, r[0] / 1000);

    // UNIT TEST 1.3 (tilt-y zp)
    n[0] = 0; n[1] = -sin(beta); n[2] = cos(beta);
    f[0] = 0.0; f[1] = 1 / T_MIN_um * 0.5;
    freq2zpXYZ(f, n, u, lambda_um, r);
    printf("UNIT TEST 1.3 (tilt-y zp): Geometric: %0.3f, fn1: %0.3f\n", rx, r[1] / 1000);

    // 1.4 UNIT TEST nominal zone plate virtual source
    beta = 0.0;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    u[0] = 0; u[1] = 0; u[2] = -p;
    fx = 1 / T_MIN_um * 0.5;
    fy = 0;
    fxy[0] = fx; fxy[1] = fy;
    freq2zpXYZ(fxy, n, u, lambda_um, r);
    rx = tan(asin(lambda_nm / T_MIN_nm / 2));
    std::cout << "UNIT TEST 1.4 (nominal zp): Geometric: " << rx << ", fn1: " << r[0] / 1000 << std::endl;

    // UNIT TEST 1.5 (tilt-x zp)
    beta = 0.1;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    freq2zpXYZ(fxy, n, u, lambda_um, r);
    printf("UNIT TEST 1.5 (tilt-x zp): Geometric: %0.3f, fn1: %0.3f\n", rx, r[0] / 1000);



    // Unit tests 2
    r[0] = tan(asin(lambda_nm / T_MIN_nm / 2)) * 1e3; 
    r[1] = 0; 
    r[2] = 1e3;
    zpCoord2KVector(r, k);



    std::cout << "UNIT TEST 2.1 (coord to freq): Assert geometric: " << k[0] << " = fn1: " << 1/T_MIN_um/2 * lambda_um << std::endl;

    // Unit tests 3
    beta = 0.1;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    u[0] = 0; u[1] = 0; u[2] = p;
    fx = 1 / T_MIN_um * 0.5;
    fy = 0;
    fxy[0] = fx; fxy[1] = fy;

    freq2zpXYZ(fxy, n, u, lambda_um, r);

    bz[0] = n[0]; bz[1] = n[1]; bz[2] = n[2];
    by[0] = 0; by[1] = 1; by[2] = 0;
    crossProduct(n, by, bx);

    zpXYZ2UxUy(r, u, bx, by, bz, U);
    std::cout << "UNIT TEST 3.1 (nominal zp): Assert bz: " << U[2] << " = 0" << std::endl;

    // Unit test 4
    zpUxUy2XYZ(U, u, bx, by, bz, r2);
    double diff = sqrt((r2[0]-r[0])*(r2[0]-r[0]) + (r2[1]-r[1])*(r2[1]-r[1]) + (r2[2]-r[2])*(r2[2]-r[2]));
    std::cout << "UNIT TEST 4.1 (nominal zp): Assert |r2 - r1| = 0: " << diff << std::endl;



    // Unit test 5
    r[0] = tan(asin(lambda_nm / T_MIN_nm / 2)) * 1e3; 
    r[1] = 0; 
    r[2] = 1e3;
    u[0] = 0; u[1] = 0; u[2] = p;


    double opl = xyz2OPL(r, u, q, lambda_um);
    double opl0 = xyz2OPL(u, u, q, lambda_um);
    std::cout << "UNIT TEST 5.1 (nominal zp): Assert OPD = 1133.048 = " << (opl - opl0) << std::endl;

    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    u[0] = 0; u[1] = 0; u[2] = p;
    fx = 1 / T_MIN_um * 0.5;
    fy = 0;
    fxy[0] = fx; fxy[1] = fy;
    freq2zpXYZ(fxy, n, u, lambda_um, r);
    rx = tan(asin(lambda_nm / T_MIN_nm / 2));

    beta = 0.1;
    n[0] = -sin(beta); n[1] = 0; n[2] = cos(beta);
    freq2zpXYZ(fxy, n, u, lambda_um, r);

    opl = xyz2OPL(r, u, q, lambda_um);
    opl0 = xyz2OPL(u, u, q, lambda_um);

    std::cout << "UNIT TEST 5.1 (nominal zp): Assert OPD = 1133.048 = " << (opl - opl0) << std::endl;




    return 0;
}