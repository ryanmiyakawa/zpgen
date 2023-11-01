/**
 * 
 *  zpGenPlus.cpp
 *  ZPGenPlus
 *  Created by Ryan Miyakawa on 2/17/17.
 *  Copyright © 2017 Ryan Miyakawa and Henry Wang. All rights reserved.
 * 
 */

#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string.h>

#include "zpUtils.h"
#include "patternFileUtils.h"

using namespace std;




double objectiveFn(double R, double th, double N, double * p, double q0, double * bx, double * by,
                   double phase, double lambda)
{
    // Compute xyz space coordinates from (r,th)
    // U are the bx,by coordinates (in-plane coordinates)
    double * U = new double[3];
    U[0] = R * cos(th);
    U[1] = R * sin(th);
    U[2] = 0;

    // Covert these back to xyz using basis vectors
    double * r = new double[3];
    
    zpUxUy2XYZ(U, p, bx, by, r);

    // Compute OPL in waves
    double opd = xyz2OPL(r, p, q0, lambda) - xyz2OPL(p, p, q0, lambda);
    double zpTerm = -N / 2;

    return opd + zpTerm + phase;
}

double secantSolve(double dRGuess, double th, double N, double * p, double q0, double * bx, double * by, double phase, double lambda)
{
    double tolX = 0.00001;
    int maxIter = 50;

    //Make guesses:
    double R1, R2, fR1, fR2, R0;
    R1 = dRGuess;
    R2 = R1 * 1.02;
    for (int currentIter = 0; currentIter < maxIter; currentIter++)
    {
        fR1 = objectiveFn(R1, th, N, p, q0, bx, by, phase, lambda);
        fR2 = objectiveFn(R2, th, N, p, q0, bx, by, phase, lambda);

        R0 = R1 - fR1 * (R1 - R2) / (fR1 - fR2);

        // Check convergence
        if (abs(R0 - R1) < tolX) 
            return R0;  // Converged!

        //Set new guesses
        R2 = R1;
        R1 = R0;
    }
    printf("MAXIMUM ITERATIONS REACHED\n");
    return R1;
}





void makeZP(
    double zTol, 
    double lambda_nm, 
    double * p, 
    double q,  
    double * k_0, 
    double * bx, 
    double * by, 
    double obscurationSigma, 
    double NA, 
    int nZerns, 
    double *orders, 
    int customMaskIdx, 
    double anamorphicFac, 
    double ZPCPhase, 
    double APD, 
    double APD_window, 
    double ZPCR2, 
    double ZPCR1, 
    double bias_nm, 
    int File_format, 
    int Opposite_Tone, 
    int FSIdx, 
    double buttressGapWidth, 
    double buttressPeriod, 
    long block_size, 
    int NoP, 
    int IoP, 
    double blockGrid_pm, 
    int layerNumber, 
    int nwaUnitSelection,
    char * fileName)
{

    clock_t begin = clock();
    long block_unit_size_pm = 500;
    double blockGrid_um = blockGrid_pm/1000000;
    double lambda, bias_um;
    lambda = lambda_nm / 1000;
    bias_um = bias_nm/1000;

    // Define plane normal
    double * n_hat = new double[3];
    crossProduct(bx, by, n_hat);

    // Need p and q to be positive
    // Need p to define NA, assume NA is always defined on fast side
    int virtualObject = 0;

    if (p[2] < 0)
    {
        virtualObject = 1;
        p[0] = abs(p[0]);
        printf("Detected Virtual Object");
    }
    

    // For WRV or GTX  block unit is passed in as nwa unit
    if (File_format == 3 || File_format == 4) {
        block_unit_size_pm = nwaUnitSelection; 

        printf("block_unit_size_pm = %ld", block_unit_size_pm);
    }
    

    
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


    double blockGrid_pixels = 0;
    int computeDoseByGriddedAreaFraction = 0;

    // WRV compute blockGrid in units of pixels:
    if (blockGrid_pm >= block_unit_size_pm && blockGrid_pm > 0){
        blockGrid_pixels = blockGrid_pm /   ((double) block_unit_size_pm);
        printf("Rounding WRV pixels by \t\t= %d pixels\n", (int)roundf(blockGrid_pixels));

        computeDoseByGriddedAreaFraction = 1;
    } else {
        blockGrid_pixels = 0;
    }
    
    
    printf("NA \t\t\t= %f \n", NA);
    printf("lambda \t\t\t= %0.3f nm\n", lambda_nm);
    printf("p \t\t\t= [%0.3f, %0.3f, %0.3f] um\n", p[0], p[1], p[2]);
    printf("image length \t\t= %g um\n", q);

    printf("ZP Normal \t\t\t= [%0.3f, %0.3f, %0.3f] \n", n_hat[0], n_hat[1], n_hat[2]);
    printf("Incident k-vector \t= [%0.3f, %0.3f, %0.3f] \n", k_0[0], k_0[1], k_0[2]);
    printf("ZP basis \t\t= {[%0.3f, %0.3f, %0.3f], [%0.3f, %0.3f, %0.3f]} \n", bx[0], bx[1], bx[2], by[0], by[1], by[2]);


    //Aberration order
    for (int j = 0; j < nZerns; j++)
    {
        printf("Zernike Number \t\t= %f \n", orders[j]);
        printf("Weight \t\t\t= %f waves\n", orders[j + nZerns]);
    }


    printf("WRV block grid \t\t = %g\n", blockGrid_pm);
    printf("BlockGrid pixels \t = %d\n", (int)roundf(blockGrid_pixels));

    
    printf("Max zone phase error\t= lambda/%d\n", (int)(1 / zTol));
    printf("ZPC phase shift \t= %f in radians\n", ZPCPhase);
    printf("Apodized term \t\t= %f %% transmission\n", APD * 100);
    if (APD_window != 0)
    {
        printf("Apodized window \t= %f \n", APD_window);
    }
    printf("Zone width bias \t= %f nm\n", bias_nm);
    printf("Is Free Standing? \t= %d \n", FSIdx);
    printf("W\t\t\t= %f dr \n", buttressGapWidth);
    printf("T\t\t\t= %f dr \n", buttressPeriod);
    printf("Block size for wrv file: %ld\n", block_size);
    printf("Generating partial %d of total %d total.\n", IoP, NoP);
    if (Opposite_Tone == 1)
    {
        printf("Tone swapped");
    }
    int wrv_split;

    printf("\n");

    if (File_format == 0)
    {
        printf("Building NWA for nanowriter\n");
    }
    else if (File_format == 1)
    {
        printf("GDS + text file \n");
    }

    else if (File_format == 2)
    {
        printf("Save as GDSII file \n");
    }
    else if (File_format == 3)
    {
        printf("Building WRV file \n");
    }
    else if (File_format == 4)
    {
        printf("Building GTX file\n");
    }
    else
    {
        printf("unknown file format %d\n", File_format);
    }

    if (NoP == 1 && IoP != 1)
    {
        printf("\nIndex of part has to be 1 if number of parts is 1 !\n");
        system("PAUSE");
    }
    if (nwaUnitSelection > 0)
    {
        printf("NWA unit selection is %d with value %0.5f.\n", nwaUnitSelection, nwaUnit);
    }

    printf("Output file: %s\n", fileName);

    /*************************/
    // BEGIN NEW ZP CODE:

    FILE *outputFile = NULL;
    FILE *supportTextFile = NULL;
    double dbscale = 0; // db unit to microns
    int numBlocksOnSide = 1;
    double rGuess, rGuessp1, Rn, Rnp1, dr, buttressWidth, alphaBT, alphaZT, alpha, 
                x, y, cx, cy, R1, R2, dR1, startAngle, currentAngle, arcStart, phase, 
                RCM = 0, tR1, tR2, f, pNA, RN, rNp, rNq, RNp1, dRN, 
                Rnpa2, Rnma2, tR1pa, tR1ma, tR2pa, tR2ma, tR1nb, tR2nb, tR1panb, tR1manb, tR2panb, tR2manb,
                minRelDose, maxRelDose, doseBias, minAreaFraction, maxAreaFraction, minClockSpeed, maxClockSpeed,  zpTiltOffsetX, zpTiltOffsetY,
                zpCenterX, zpCenterY, offsetX = 0, offsetY = 0, drawAngle, ux, uy;

    double * r = new double[3]; 
    double * u = new double[3];
    double * PC = new double[2];

    long totalPoly = 0;
    unsigned char gdsPost[8];
    unsigned char polyPre[16];
    unsigned char polyPost[4];
    unsigned char polyForm[4];

    // =============================================================//

    /**
     * @brief Compute max zone number
     * Use na and lambda to compute d_min.  Then generate a circle in k-space around k0, and compute the max zone number from these points
     */

    double T_min = lambda / NA;
    int N_max = 1;
    int N_min = 1;
    double th = 0;
    double fx = 0, fy = 0;
    double opd = 0;
    int tempN_min;
    int tempN_max;
    double * fq = new double[2];
    int numPoints = 200;
    for (int i = 0; i < numPoints; i++){
        th = i * 2 * M_PI / numPoints;
        fq[0] = 1 / T_min * cos(th) + k_0[0] / lambda;// is this right?
        fq[1] = 1 / T_min * sin(th) + k_0[1] / lambda;

        // find opd for each:
        freq2zpXYZ(fq, n_hat, p, lambda, r);
        opd = (double)((long double)xyz2OPL(r, p, q, lambda) - (long double)xyz2OPL(p, p, q, lambda));

        // printf("th: %0.2f, opd: %0.3f\n", th, opd);
        
        tempN_min = (int) 2 * opd;
        tempN_max = (int) 2*opd + 1;

        if (tempN_min < N_min){
            N_min = tempN_min;
        }
        if (tempN_max > N_max){
            N_max = tempN_max;
        }
       
    }

    printf("N_max: %d, N_min: %d\n", N_max, N_min);

    if (File_format == 3 || File_format == 4)
    {                                  // WRV: shift by one half block size
        numBlocksOnSide = layerNumber; // For WRV files we send num blocks on the GDSLayer input

        offsetX += numBlocksOnSide * block_size * block_unit_size_pm / 2 / 1000000;
        offsetY += numBlocksOnSide * block_size * block_unit_size_pm / 2 / 1000000;
    }

    if (File_format == 4)
    { // GTX: shift by one half block size
        offsetX +=  block_size / 2 ;
        offsetY +=  block_size / 2 ;
    }

    // Compute dRN and dR1 to bound dose information
    f          = 1/(1/norm2(p) * std::copysign(1.0, p[2]) + 1/q);
    RN          = secantSolve(sqrt(N_max*lambda*f), 0, N_max, p, q, bx, by, 0, lambda);
    RNp1        = secantSolve(sqrt((N_max + 1)*lambda*f), 0, N_max + 1, p, q, bx, by, 0, lambda);
    dRN         = RNp1 - RN;
    R1          = secantSolve(sqrt(N_min*lambda*f), 0, N_min, p, q, bx, by, 0, lambda);
    R2          = secantSolve(sqrt((N_min + 1)*lambda*f), 0, N_min + 1, p, q, bx, by, 0, lambda);
    dR1         = R2 - R1;

    // Define area fraction (nominal area)/(biased or new area)
    if (computeDoseByGriddedAreaFraction){ 
        // Smallest area fraction
        minAreaFraction = 1;

        // Largest possible area fraction
        maxAreaFraction  = (dRN + blockGrid_um)/(dRN - bias_um - blockGrid_um);
    } else {
        maxAreaFraction     = dRN / (dRN - bias_um); // furthest from 1 to 1 i.e. (1.24)
        minAreaFraction     = dR1 / (dR1 - bias_um); // closest to 1 i.e. (1.02)
    }

    maxRelDose = maxAreaFraction;
    minRelDose = minAreaFraction;

    
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
            
            initWRV(outputFile, minRelDose, maxRelDose, block_size, block_unit_size_pm);
            break;
        case 4: // GTX
            dbscale = 1000000/block_unit_size_pm; // db unit to microns
            if ((outputFile = fopen(strcat(fileName, ".gtx"), "wb")) == NULL)
                printf("Cannot open file.\n");
            
            initGTX(outputFile, block_unit_size_pm/1000000, (long)offsetX, (long)offsetY);
            break;
    }

    long clockSpeed;
    bool isGapZone = false;

    // If we're printing obscuration only, then set N to 0:
    int N;
    if (customMaskIdx == 16)
        N = 0;

    // Initialze guesses
    rGuess = sqrt(N_min*lambda*f);
    rGuessp1 = sqrt((N_min + 1)*lambda*f);

    for (int n = N_min; n <= N_max; n++)
    {
        if (n == N_max){
            printf("max n");
        }
        // Compute initial R at an arbitrary angle.  This will seed information regarding RoC for zone tol:
        Rn          = secantSolve(rGuess, 0, n, p, q, bx, by, 0, lambda);
        Rnp1        = secantSolve(rGuessp1, 0, n + 1, p, q, bx, by, 0, lambda);
        dr         = Rnp1 - Rn;

        // Use current zone radii for new guesses
        rGuess     = rGuessp1;
        rGuessp1   = rGuessp1 * sqrt(((double)n + 1)/((double)n));
        
    
        // If buttresses are specified, condition segment width on zone parity
        buttressWidth = 0;
        // ButtressWidth is the width of the buttressed zone segment. Gap is the space in between - technically, the actual buttress
        switch (FSIdx)
        {
        case 0: // % no buttresses
            if (n % 2 == Opposite_Tone) { continue; }
            buttressWidth = 0;
            break;
        case 1: // % gapped zones
            if (n % 2 == Opposite_Tone) { continue;}
            buttressWidth = dr * (buttressPeriod - buttressGapWidth);
            break;
        case 2: // full zones + gaps
            if (n % 2 == Opposite_Tone)
            {               //% Even zones get gaps
                if (n == N) // Don't make last zone of gaps
                    continue;

                isGapZone = true;
                buttressWidth = dr * (buttressGapWidth);
            }
            else
            { //% Odd zones get Full zones
                isGapZone = false;
                buttressWidth = 0;
            }
            break;
        }

        /* We can comply with zone tol constraint by leveraging the radius of curvature of the present zone programmatically. */
        alphaZT = 2 * acos(1 / (zTol * lambda / Rn + 1)); // Full angle subtended by valid arc segment

        // Compare with full angle specified by buttress width
        alphaBT = buttressWidth / Rn;

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
        startAngle = ((double)(rand() % 1000)) / 1000 * 2 * M_PI; // randomize start
        currentAngle = startAngle;

        double gapZoneSize = 0;


        // Loop through angles
        while (true){

            if (currentAngle >= startAngle + 2 * M_PI)
                break;
        /**
         * @brief 
         * 
         *
         * 1) Using current angle, compute baseline Rn for n (no added phase component)
         * 2) Use computed Rn and currentAngle to compute [cx,cy] to figure out phase and amplitude
         * 3) Recompute Rn_p and Rnp1_p with new phase and [cx,cy] with new phase for pupil mask
         */

            
        //  1) Using current angle, compute baseline Rn for n (no added phase component)
            Rn          = secantSolve(Rn, currentAngle, n, p, q, bx, by, 0, lambda);
            Rnp1        = secantSolve(Rnp1, currentAngle, n + 1, p, q, bx, by, 0, lambda);

        //  2) Use computed Rn and currentAngle to compute [cx,cy] to figure out phase and amplitude
            zpRTh2PCCxCy(Rn, currentAngle, k_0, p, bx, by, lambda, NA, PC);
            cx = PC[0];
            cy = PC[1]/anamorphicFac;

            // Compute phase terms in waves from Zernike polynomials and custom phase
            phase = getPhaseTerm(cx, cy, orders, nZerns, ZPCPhase, ZPCR1, ZPCR2);
            phase += customPhase(cx, cy, customMaskIdx) / (2 * M_PI);

            // Accept or reject trap based on pupil boundaries, obscuration, and custom mask     
            if (!bIsInGeometry(cx, cy, obscurationSigma) || 
            (customMaskIdx != 0 && !bIsInCustomMask(cx, cy, customMaskIdx)))
            {
                currentAngle = currentAngle + alpha;
                continue;
            }

        // * 3) Recompute Rn_p and Rnp1_p and [cx, cy] with new phase
            Rn          = secantSolve(Rn, currentAngle, n, p, q, bx, by, phase, lambda);
            Rnp1        = secantSolve(Rnp1, currentAngle, n + 1, p, q, bx, by, phase, lambda);

            // Compute CM to optimized trap to arc and trap coords
            dr = Rnp1 - Rn;
            RCM = (Rnp1 * Rnp1 * Rnp1 - Rn * Rn * Rn) / (Rnp1 * Rnp1 - Rn * Rn) * 2 / 3 * sin(alpha) / alpha;
            //RCM = (Rn + dr/2)*sin(alpha)/alpha; // CM of arc: center trap on arc CM rather than matching

            Rnpa2          = secantSolve(Rn, currentAngle + alpha/2, n, p, q, bx, by, phase, lambda);
            Rnma2          = secantSolve(Rn, currentAngle - alpha/2, n, p, q, bx, by, phase, lambda);


            drawAngle = currentAngle;
            if (File_format == 4){ // GTX
                doseBias = dr / (dr - bias_um); // Inverse of zone area bias
                if (isGapZone)
                {
                    exportArcGTX(Rn - bias_um / 2, dr + bias_um, drawAngle, alpha, doseBias, ((double)(block_unit_size_pm))/1000000, 0, 0, outputFile);
                }
                else
                {
                    exportArcGTX(Rn + bias_um / 2, dr - bias_um, drawAngle, alpha, doseBias, ((double)(block_unit_size_pm))/1000000, 0, 0, outputFile);
                }
            }
            else { // GDS or WRV format
                
                // Unbiased CM radii
                tR1nb = (RCM - dr/2)*1/cos(alpha/2);
                tR2nb = (RCM + dr/2)*1/cos(alpha/2);

                // Appy biases
                if (isGapZone){
                    tR1 = tR1nb - bias_um/2;
                    tR2 = tR2nb + bias_um/2;
                } else {
                    tR1 = tR1nb + bias_um/2;
                    tR2 = tR2nb - bias_um/2;
                }
                
                // Compute corrected radii at the edges of trap instead of center:
                tR1ma = tR1 + (Rnma2 - Rn);
                tR1pa = tR1 + (Rnpa2 - Rn);
                tR2ma = tR2 + (Rnma2 - Rn);
                tR2pa = tR2 + (Rnpa2 - Rn);

                // version with no bias
                tR1manb = tR1nb + (Rnma2 - Rn);
                tR1panb = tR1nb + (Rnpa2 - Rn);
                tR2manb = tR2nb + (Rnma2 - Rn);
                tR2panb = tR2nb + (Rnpa2 - Rn);
  
                // Coordinates of trap in real distance
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
                double trapCoords[10];
                    trapCoords[0] =  dbscale*trapCoords_um[0];
                    trapCoords[1] =  dbscale*trapCoords_um[1];
                    trapCoords[2] =  dbscale*trapCoords_um[2];
                    trapCoords[3] =  dbscale*trapCoords_um[3];
                    trapCoords[4] =  dbscale*trapCoords_um[4];
                    trapCoords[5] =  dbscale*trapCoords_um[5];
                    trapCoords[6] =  dbscale*trapCoords_um[6];
                    trapCoords[7] =  dbscale*trapCoords_um[7];
                    trapCoords[8] =  dbscale*trapCoords_um[0];
                    trapCoords[9] =  dbscale*trapCoords_um[1];


                if (computeDoseByGriddedAreaFraction){
                    // compute clockspeed by area fraction, (not working right now):
                    clockSpeed = getPolyClockSpeedFromRadiusRatio(minAreaFraction, maxAreaFraction, dr, trapCoords, dbscale, blockGrid_pixels);

                } else {
                    // Old way to compute clock speed, still working
                    clockSpeed = getPolyClockSpeed(dR1, dRN, dr, bias_um);
                 }
            
                // Export shape
                exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format, clockSpeed, blockGrid_pixels);
            }

            trapCount++;

            // Increment angle by amount depending on whether we require a gap or another trap to satisfy zTol requirement
            if (buttressWidth == 0)
            {
                currentAngle = currentAngle + alpha;
            }
            else if (isGapZone)
            {
                // If the bridges in the odd zones require more than one trap to make, then keep a running tab of size
                gapZoneSize += alpha;
                if (gapZoneSize >= alphaBT)
                {
                    currentAngle = currentAngle + dr * buttressPeriod / RCM - alphaBT;
                    gapZoneSize = 0;
                }
                else
                {
                    currentAngle = currentAngle + alpha;
                }
            }
            else
            {
                if (alphaBT > alphaZT)
                {
                    if (arcStart < 0)
                    { // This segment requires multipe traps to satsify ztol.  Start counting
                        arcStart = currentAngle;
                    }
                    if ((currentAngle - arcStart + alpha) > alphaBT)
                    { // Require a gap here, reset counter
                        currentAngle = currentAngle + alpha + dr * buttressGapWidth / RCM;
                        arcStart = -1;
                    }
                    else
                    { // keep adding traps in this segment
                        currentAngle = currentAngle + alpha;
                    }
                }
                else
                {
                    currentAngle = currentAngle + dr * buttressPeriod / RCM; // This segment satisfies ztol
                }
            }

        } // End angle loop

        totalPoly += trapCount;
        if (File_format == 0)
            printf("Finished zone %d with %ld arcs\n", n, trapCount);
        else
            printf("Finished zone %d with %ld traps.  \tR_%d = %0.5f \n", n, trapCount, n, Rn);


    } // End zone loop

    // Write obscuration.  Note: cannot write obscuration in arcs presently
    if (Opposite_Tone && obscurationSigma > 0 && File_format != 0 && File_format != 4)
    {
        printf("Writing obscuration with sigma: %0.2f\n", obscurationSigma);
        printf("Obscuration anamorphic factor: %0.2f\n", anamorphicFac);
        if (anamorphicFac != 1)
        {
            printf("Writing asymmetric obscuration\n");
        }

        int nObsPts = 360;
        double qXCoords[nObsPts];
        double qYCoords[nObsPts];

        double thStep = 2 * M_PI / nObsPts;
        double theta;
        for (int i = 0; i < nObsPts; i++)
        {
            theta = thStep * i + M_PI / 2;

            // Convert pupil coordinates to zone plate coordinates
            fq[0] = (NA * obscurationSigma * cos(theta) + k_0[0]) / lambda;
            fq[1] =( NA * obscurationSigma * sin(theta) + k_0[1]) / lambda;

            freq2zpUxUy(fq, n_hat, p, bx, by, lambda, u);

            qXCoords[i] = u[0];
            qYCoords[i] = u[1];
        }

        double azRot = 0;
        for (int k = 0; k < (nObsPts / 2 - 1); k++)
        {
            totalPoly++;

            double trapCoords[10];
            trapCoords[0] =  dbscale*(offsetX + cos(azRot)*qXCoords[k] - sin(azRot) * qYCoords[k]);
            trapCoords[1] =  dbscale*(offsetY + sin(azRot)*qXCoords[k] + cos(azRot) * qYCoords[k]);

            trapCoords[2] =  dbscale*(offsetX + cos(azRot)*qXCoords[k+1] - sin(azRot) * qYCoords[k+1]);
            trapCoords[3] =  dbscale*(offsetY + sin(azRot)*qXCoords[k+1] + cos(azRot) * qYCoords[k+1]);

            trapCoords[4] =  dbscale*(offsetX + cos(azRot)*qXCoords[nObsPts - k - 2] - sin(azRot) * qYCoords[nObsPts - k - 2]);
            trapCoords[5] =  dbscale*(offsetY + sin(azRot)*qXCoords[nObsPts - k - 2] + cos(azRot) * qYCoords[nObsPts - k - 2]);

            trapCoords[6] =  dbscale*(offsetX + cos(azRot)*qXCoords[nObsPts - k - 1] - sin(azRot) * qYCoords[nObsPts - k - 1]);
            trapCoords[7] =  dbscale*(offsetY + sin(azRot)*qXCoords[nObsPts - k - 1] + cos(azRot) * qYCoords[nObsPts - k - 1]);

            trapCoords[8] =  dbscale*(offsetX + cos(azRot)*qXCoords[k] - sin(azRot) * qYCoords[k]);
            trapCoords[9] =  dbscale*(offsetY + sin(azRot)*qXCoords[k] + cos(azRot) * qYCoords[k]);
        
            // Export shape
            exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format, 65535, blockGrid_pixels);

        }
    }

    switch (File_format)
    {
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
    case 4: // GTX
        renderGTX(outputFile);
        break;
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Finished zone plate with %ld shapes\n", totalPoly);
    printf("\nZP Generation took: %0.3f seconds\n", elapsed_secs);


}


// Primary Gateway
int main(int argc, char **argv){
    // Geometric vectors:
    double * p = new double[3]; // Vector from origin to zone plate along optical axis
    double * n = new double[3]; // Unit vector normal to zone plate
    double * k_0 = new double[3]; // Unit vector in propagation direction
    double * bx = new double[3]; // Zone plate basis vector x
    double * by = new double[3]; // Zone plate basis vector y

    double zTol, lambda_nm, q, obscurationSigma, NA, anamorphicFac, ZPCPhase, APD, APD_window, ZPCR2, ZPCR1, bias_nm,
    buttressGapWidth, buttressPeriod, blockGrid_pm;
    int nZerns, customMaskIdx, File_format, Opposite_Tone, FSIdx, layerNumber, nwaUnitSelection, NoP, IoP;
    long block_size;
    char * fileName;




    if (argc == 1)
    {
        printf("No input parameters!!! \n");
        system("PAUSE");

        double zeroTriple[] = {0, 0, 0};
        double zeroDouble[] = {0, 0};
        double zeroEmpty[] = {};
        double tempP[] = {0, 0, 100};


        double tempBxFlat[] = {1, 0, 0};
        double tempBxTilt[] = {cos(6 * M_PI/180), 0, -sin(6 * M_PI/180)};
        double tempByFlat[] = {0, 1, 0};
        double tempByTilt[] = { 0, cos(6 * M_PI/180), -sin(6 * M_PI/180)};


        double tempK[] = {0, 0, 1};
        double tempKTilt[] = {sin(6 * M_PI/180), 0, cos(6 * M_PI/180)};

        double * bx;
        double * by;
        char str1[] = "debug0";
        char str2[] = "debugTilt";
        char str3[] = "debugOAxis";

        char * fileName;

        int setting = 3;
        if (setting == 1){
            bx = tempBxFlat;
            by = tempByFlat;
            k_0 = tempK;
            fileName = str1;
        } else if (setting == 2) {
            bx = tempBxFlat;
            by = tempByTilt;
            fileName = str2;
            k_0 = tempK;
        } else if (setting == 3) {
            // CRA
            k_0 = tempKTilt;
            bx = tempBxFlat;
            by = tempByFlat;
            fileName = str3;
        }

        
        makeZP(
                0.01, //zTol, 
                13.5, //lambda_nm, 
                tempP, // p
                100000, // q
                k_0, //k, 
                bx, //bx, 
                by, //by, 
                0, //obscurationSigma, 
                0.02, //NA, 
                0, //nZerns, 
                zeroEmpty, //orders, 
                0, //customMaskIdx, 
                1, //anamorphicFac, 
                0, //ZPCPhase, 
                0, //APD, 
                0, //APD_window, 
                0, //ZPCR2, 
                0, //ZPCR1, 
                10, //bias_nm, 
                1, // GDS File_format, 
                0, // Opposite_Tone, 
                0, //FSIdx, 
                0, //buttressGapWidth, 
                0, //buttressPeriod, 
                0, //block_size, 
                0, //NoP, 
                0, //IoP, 
                0, //blockGrid_pm, 
                0, //layerNumber, 
                0, //nwaUnitSelection, 
                fileName //fileName
        );


        return 0;
    }

printf("Reading input params\n");


    // return 0;


    char **argv_test;
    argv_test = argv + 1;


    zTol             = atof(*(argv_test++));
    lambda_nm        = atof(*(argv_test++));

    p[0]                   = atof(*(argv_test++));
    p[1]                   = atof(*(argv_test++));
    p[2]                   = atof(*(argv_test++));
    q                       = atof(*(argv_test++));

    k_0[0]                   = atof(*(argv_test++));
    k_0[1]                   = atof(*(argv_test++));
    k_0[2]                   = atof(*(argv_test++));
    bx[0]                  = atof(*(argv_test++));
    bx[1]                  = atof(*(argv_test++));
    bx[2]                  = atof(*(argv_test++));
    by[0]                  = atof(*(argv_test++));
    by[1]                  = atof(*(argv_test++));
    by[2]                  = atof(*(argv_test++));


    obscurationSigma = atof(*(argv_test++));
    NA               = atof(*(argv_test++)); 
    nZerns              = atoi(*(argv_test++));

    double orders[2 * nZerns];
    for (int i = 0; i < nZerns; i++)
    {
        orders[i] = atof(*(argv_test++));
        orders[i + nZerns] = atof(*(argv_test++));
    }

    customMaskIdx       = atoi(*(argv_test++)); 
    anamorphicFac       = atof(*(argv_test++));
    ZPCPhase            = atof(*(argv_test++));
    APD                 = atof(*(argv_test++));
    APD_window          = atof(*(argv_test++));
    ZPCR2               = atof(*(argv_test++));
    ZPCR1               = atof(*(argv_test++));
    bias_nm             = atof(*(argv_test++));
    File_format         = atof(*(argv_test++));
    Opposite_Tone       = atof(*(argv_test++));
    FSIdx               = atof(*(argv_test++));
    buttressGapWidth    = atof(*(argv_test++));
    buttressPeriod      = atof(*(argv_test++));
    //1 million pixel with a pixel size 500 pm (0.5 nm)
    // Thus it transfer to 500 um block size
    block_size          = atol(*(argv_test++));
    NoP                 = atoi(*(argv_test++));
    IoP                 = atoi(*(argv_test++));
    blockGrid_pm        = atof(*(argv_test++)); // block grid (WRV)


    layerNumber         = atoi(*(argv_test++));
    nwaUnitSelection    = atoi(*(argv_test++));
    fileName            = *argv_test;


    makeZP(zTol, 
        lambda_nm, 
        p, 
        q, 
        k_0, 
        bx, 
        by, 
        obscurationSigma, 
        NA, 
        nZerns, 
        orders, 
        customMaskIdx, 
        anamorphicFac, 
        ZPCPhase, 
        APD, 
        APD_window, 
        ZPCR2, 
        ZPCR1, 
        bias_nm, 
        File_format, 
        Opposite_Tone, 
        FSIdx, 
        buttressGapWidth, 
        buttressPeriod, 
        block_size, 
        NoP, 
        IoP, 
        blockGrid_pm, 
        layerNumber, 
        nwaUnitSelection, 
        fileName);

    return 0;
}
