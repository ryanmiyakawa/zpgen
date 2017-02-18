//
//  zpGenPlus.cpp
//  ZPGenPlus
//
//  Created by Ryan Miyakawa on 2/17/17.
//  Copyright © 2017 Ryan Miyakawa and Henry Wang. All rights reserved.
//
// Refactoring ZPGen for speed and accuracy improvements.  Goal is to make zone plate GDS file generation take less than 5 seconds each

#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

// Recursively computes the binomial coefficient nCk:
double nChooseK(double n, double k) {
    if (k == 0 && n >= 0)
        return 1;
    else if (n == 0 && k > 0)
        return 0;
    else
        return nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
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
        
        Rv += kParity*nChooseK((double)(n - k), (double)k)
        *nChooseK((double)(n - 2 * k), (double)(p - k))
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

void exportPolygon(long * coords, unsigned char * polyPre, unsigned char * polyPost, unsigned char * polyForm, FILE * outputFile, int File_format){
    switch (File_format){
        case 0:
            break;
        case 1:
            break;
        case 2:
            char cCoords[40];
            fwrite(polyPre, sizeof(char), 16, outputFile);
            fwrite(polyForm, sizeof(char), 4, outputFile);
            encodePoly32(coords, cCoords);
            fwrite(cCoords, sizeof(char), 40, outputFile);
            fwrite(polyPost, sizeof(char), 4, outputFile);
            break;
        case 3:
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

bool bIsInGeometry(double cx, double cy, double obscurationSigma){
    double r = sqrt(cx*cx + cy*cy);
    return (r >= obscurationSigma && r <= 1.002);
}

bool bIsInCustomMask(double cx, double cy, int customMaskIdx){
    double dx, dy;
    switch (customMaskIdx){
        case 1: // Intel MET
            return (sqrt(square(cx - .65*cos(7.5*M_PI/180)) + square(cy - .65*sin(7.5*M_PI/180))) <= .35) ||
                        (sqrt(square(cx - .65*cos(127.5*M_PI/180)) + square(cy - .65*sin(127.5*M_PI/180))) <= .35) ||
                        (sqrt(square(cx - .65*cos(247.5*M_PI/180)) + square(cy - .65*sin(247.5*M_PI/180))) <= .35);
            break;
            
        case 2: // TDS ZP 2
            dx = cx * 0.1515;
            dy = cy * 0.1515;
            return square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            break;
            
        case 3: // TDS ZP 3
            dx = cx * 0.2396;
            dy = cy * 0.2396;
            return square(dx)/square(0.2396) + square(dy)/square(0.1515) <= 1 &&  // outer ellipse
                    square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            break;
            
        case 4: // TDS ZP 4
            dx = cx * 0.285;
            dy = cy * 0.285;
            return square(dx) + square(dy + 0.1515) <= square(0.303) &&  // outer parent
                    dy + 0.1515 >= 0 &&  // x-axis
                    (dy + 0.1515) >= tan(37*M_PI/180)*(dx - 0.138) && // diagonal lines
                    (dy + 0.1515) >= -tan(37*M_PI/180)*(dx + 0.138) &&
                    square(dx) + square(dy + 0.1515) > square(0.018); // inner circle
            
            break;
        case 5: // Square aperture
            
            break;
    }
    return true;
}

double objectiveFn(double r, double th, double N, double p, double q,
                   double phase, double lambda, double beta) {
    
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
    
    return plTermP + plTermQ + zpTerm + phTerm;
}

double secantSolve(double dRGuess, double th, double N, double p, double q, double phase, double lambda, double beta){
    double tolX = 0.00001;
    int maxIter = 20;
    
    //Make guesses:
    double x1, x2, fxm1, fxm2, R0;
    x1 = dRGuess;
    x2 = x1*1.02;
    for (int currentIter = 0; currentIter < maxIter; currentIter++){
        fxm1 = objectiveFn(x1, th, N, p, q, phase, lambda, beta);
        fxm2 = objectiveFn(x2, th, N, p, q, phase, lambda, beta);

        R0 = x1 - fxm1*(x1 - x2) / (fxm1 - fxm2);
        if (abs(R0 - x1) < tolX)
            return R0;
    
        //Set new guesses
        x2 = x1;
        x1 = R0;
    }
    printf("MAXIMUM ITERATIONS REACHED\n");
    return x1;
}

void initGDS(FILE * outputFile, unsigned char * gdsPost, unsigned char* polyPre, unsigned char * polyPost, unsigned char * polyForm){
    int gdspreamble[102] = { 0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
        109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
        243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
        0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
        0, 10, 6, 6, 110, 111, 110, 97, 109, 101};
    
    int gdspostamble[8]     = {0, 4, 7, 0, 0, 4, 4, 0};
    int polypreamble[16]    = {0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0};
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
    
    // Need p to define NA, assume NA is always defined on fast side
    if (p > q){ // swap
        double tempP = p;
        p = q;
        q = tempP;
    }
    
    double obscurationSigma = atof(*(argv_test++));
    double NA_Pdep = atof(*(argv_test++)); // remove this
    int nZerns = atoi(*(argv_test++));
    
    double orders[2 * nZerns];
    for (int k = 0; k < nZerns; k++) {
        orders[k] = atof(*(argv_test++));
        orders[k + nZerns] = atof(*(argv_test++));
    }
    int customMaskIdx = atoi(*(argv_test++)); // Changed alpha to custom mask
    double beta = atof(*(argv_test++));
    double CRAAz = atof(*(argv_test++))* M_PI/180;
    double CRA = atof(*(argv_test++)) * M_PI/180;
    double NA = atof(*(argv_test++));
    double ZPCPhase = atof(*(argv_test++));
    double APD = atof(*(argv_test++));
    double APD_window = atof(*(argv_test++));
    double ZPCR2 = atof(*(argv_test++));
    double ZPCR1 = atof(*(argv_test++));
    double bias_nm = atof(*(argv_test++));
    int File_format = atof(*(argv_test++));
    int Opposite_Tone = atof(*(argv_test++));
    int FSIdx = atof(*(argv_test++));
    double buttressGapWidth = atof(*(argv_test++));
    double buttressPeriod = atof(*(argv_test++));
    double TrapDoesPower = atof(*(argv_test++));
    //1 million pixel with a pixel size 500 pm (0.5 nm)
    // Thus it transfer to 500 um block size
    int block_size = atof(*(argv_test++));
    double NoP = atol(*(argv_test++));
    double IoP = atol(*(argv_test++));
    int zpID = atoi(*(argv_test++));
    int curl_on = atoi(*(argv_test++));
    char * fileName = *argv_test;
    
    double lambda, bias_um;
    double NA_P = NA + sin(CRA);
    double D =  2 * (p*q / (p + q))*tan(asin(NA_P));
    
    lambda = lambda_nm / 1000;
    bias_um = bias_nm/1000;
    
    //From block_size x resolution limit (1000pm = 1nm)
    block_size = block_size * 10;
    
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
    
    printf("titlted angle beta \t= %g in radian\n", beta);
    printf("Max zone phase error\t= lambda/%d\n", (int)(1 / zTol));
    printf("Zone plate diameter \t= %0.2f um\n", D);
    printf("Extra Phase shift \t= %f in Radian\n", ZPCPhase);
    printf("Apodized term \t\t= %f %% transmission\n", APD * 100);
    if (APD_window != 0) {
        printf("Apodized window \t= %f \n", APD_window);
    }
    printf("Bias term \t\t= %f nm\n", bias_nm);
    printf("Free Standing \t\t= %d \n", FSIdx);
    printf("W\t\t\t= %f dr \n", buttressGapWidth);
    printf("T\t\t\t= %f dr \n", buttressPeriod);
    printf("Trap Dose Power \t= %f\n", TrapDoesPower);
    printf("Block sieze for wrv file: %d\n", block_size);
    printf("Generate %f part out of total %f parts.\n", IoP, NoP);
    if (Opposite_Tone == 1) {
        printf("Print in opposie tone!");
    }
    int wrv_split;
    if (D >= 400 && File_format == 3 && CRA == 0) {
        wrv_split = 1;
        printf("On-axis Zoneplate is larger than WRV block size, split into 9 different blocks/files");
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
    
    /*************************/
    // BEGIN NEW ZP CODE:
    
    FILE * outputFile;
    double dbscale = 10000; // db unit to microns
    double rGuess, rGuessp1, Rn, Rnp1, dr, buttressWidth, alphaBT, alphaZT, alpha, x, y, cx, cy,
                                startAngle, currentAngle, arcStart, phase, RCM, tR1, tR2;
    unsigned char gdsPost[8];
    unsigned char polyPre[16];
    unsigned char polyPost[4];
    unsigned char polyForm[4];
    if ((outputFile = fopen(fileName, "wb")) == NULL)
        printf("Cannot open file.\n");

    switch (File_format){
        case 0: // arc
            
            break;
        case 1: // GDS
            initGDS(outputFile, gdsPost, polyPre, polyPost, polyForm);
            break;
        case 2://GDS + txt
            initGDS(outputFile, gdsPost, polyPre, polyPost, polyForm);
            break;
            
            break;
        case 3: // WRV
            
            break;
    }
    
    /* Figure out number of zones by counting the number of waves that separate rN and f */
    double f = 1/(1/p + 1/q);
    double pNA = NA + sin(CRA);
    double RN = p * tan(asin(pNA));    // R in plane of ZP
    double rN = sqrt(RN*RN + p*p);     // r is hypotenuse of RN and f
    
    //Compute number of zones in this zone plate using RN PLD
    int N = (int)(2*(rN - p)/lambda);

    for (int n = 1; n <= N; n++){
        //Compute initial R at an arbitrary angle.  This will seed information regarding RoC for zone tol:
        rGuess     = sqrt(n*lambda*f);
        rGuessp1   = sqrt((n + 1)*lambda*f);
        Rn         = secantSolve(rGuess, 0, n, p, q, 0, lambda, beta);
        Rnp1       = secantSolve(rGuessp1, 0, n + 1, p, q, 0, lambda, beta);
        dr         = Rnp1 - Rn;
    
    //If buttresses are specified, condition segment width on zone parity
        buttressWidth = 0;
        switch (FSIdx){
            case 0:// % no buttresses
                if (n % 2 == 0){
                    continue;
                }
                buttressWidth = 0;
                break;
            case 1:// % gapped zones
                if (n % 2 == 0){
                    continue;
                }
                buttressWidth = dr*(buttressPeriod - buttressGapWidth);
                break;
            case 2: // full zones + gaps
                if (n % 2 == 0){ //% Even zones get gaps
                    buttressWidth = dr*(buttressGapWidth);
                }else{ //% Odd zones get Full zones
                    buttressWidth = 0;
                }
                break;
        }
    
    /* We can comply with zone tol constraint by leveraging the radius of curvature of the present zone programmatically. */
    alphaZT = 2*acos(1/(zTol * lambda / Rn + 1)); // Full angle subtended by valid arc segment
    
    // Compare with full angle specified by buttress width
    alphaBT = buttressWidth/Rn;
    
    // Set alpha to the tighter of the two constraints
    if (buttressWidth != 0 && alphaBT < alphaZT)
        alpha = alphaBT;
    else
        alpha = alphaZT;
    
    //For dealing with mixed conditions above
    arcStart = -1;
    
    // Log trap count:
    long trapCount = 0;
    
    //Loop through angle
    startAngle     = ((double) (rand() % 1000))/1000 * 2 * M_PI ; // randomize start
    currentAngle = startAngle;
    
    while(true){
        if (currentAngle >= startAngle + 2*M_PI){
            break;
        }
        // Get relative coordinates
        x = Rn*cos(currentAngle);
        y = Rn*sin(currentAngle);
    
        // Convert coordinates to pupil coordinates
        if (CRA == 0){
            cx = sin(atan(x/p))/NA;
            cy = sin(atan(y/p))/NA;
        }
        else{
            cx = (sin(atan(x/p)) - sin(CRA)*cos(CRAAz))/NA;
            cy = (sin(atan(y/p)) - sin(CRA)*sin(CRAAz))/NA;
        }
    
        // Accept or reject trap based on non-aberrated geometry
        if (!bIsInGeometry(cx, cy, obscurationSigma)){
            currentAngle = currentAngle + alpha;
            continue;
        }
        
        // Apply custom mask
        if (customMaskIdx != 0 && !bIsInCustomMask(cx, cy, customMaskIdx)){
            currentAngle = currentAngle + alpha;
            continue;
        }
    
        // Compute phase terms
        phase = getPhaseTerm(cx, cy, orders, nZerns, ZPCPhase, ZPCR1, ZPCR2);
    
        // Use initial R as seeds for each subsequent comp
        Rn   = secantSolve(Rn, currentAngle, n, p, q, phase, lambda, beta);
        Rnp1 = secantSolve(Rnp1, currentAngle, n + 1, p, q, phase, lambda, beta);
        
        // Compute CM and trap coords
        dr = Rnp1 - Rn;
        RCM = (Rn + dr/2)*sin(alpha)/alpha; // CM of arc: center trap on arc CM rather than matching vertices
        tR1 = (RCM - dr/2)*1/cos(alpha/2) + bias_um/2;
        tR2 = (RCM + dr/2)*1/cos(alpha/2) - bias_um/2;
    
        long trapCoords[10];
        trapCoords[0] = (long) dbscale*tR1*cos(currentAngle - alpha/2);
        trapCoords[1] = (long) dbscale*tR1*sin(currentAngle - alpha/2);
        trapCoords[2] = (long) dbscale*tR1*cos(currentAngle + alpha/2);
        trapCoords[3] = (long) dbscale*tR1*sin(currentAngle + alpha/2);
        trapCoords[4] = (long) dbscale*tR2*cos(currentAngle + alpha/2);
        trapCoords[5] = (long) dbscale*tR2*sin(currentAngle + alpha/2);
        trapCoords[6] = (long) dbscale*tR2*cos(currentAngle - alpha/2);
        trapCoords[7] = (long) dbscale*tR2*sin(currentAngle - alpha/2);
        trapCoords[8] = (long) dbscale*tR1*cos(currentAngle - alpha/2);
        trapCoords[9] = (long) dbscale*tR1*sin(currentAngle - alpha/2);
    
        // Export shape
        exportPolygon(trapCoords, polyPre, polyPost, polyForm, outputFile, File_format);
    
        trapCount++;
    
        // Increment angle by amount depending on whether we require a gap or another trap to satisfy zTol requirement
        if (buttressWidth == 0){
            currentAngle = currentAngle + alpha;
        }
        else{
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
        
        }
        printf("Finished zone %d with %ld traps\n", n, trapCount);
    }
    
    switch (File_format){
        case 0: // arc
            
            break;
        case 1: // GDS
            renderGDS(outputFile, gdsPost);
            break;
        case 2: // GDS + txt
            
            
            break;
        case 3: // WRV
            
            break;
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("\nZP Generation took: %0.3f seconds\n", elapsed_secs);
    system("PAUSE");
    
    
    return 0;
}




