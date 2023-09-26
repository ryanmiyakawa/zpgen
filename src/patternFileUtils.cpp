

#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string.h>

#include "zpUtils.h"

using namespace std;


// Computes dose ratio from zone width bias
double getDoseFromBias(double dr, double bias)
{
    return dr / (dr - bias);
}

void writeSupportTxtHeader(double dR1, double dRN, double bias, FILE *supportTextFile)
{
    fprintf(supportTextFile, "Dose %f %f\n", getDoseFromBias(dRN, bias), getDoseFromBias(dR1, bias));
}

void writeSupportTxtZoneDose(long clockSpeed, int n, double Rn, FILE *supportTextFile)
{
    fprintf(supportTextFile, "Zone %d %ld %f\n", n, clockSpeed, Rn);
}

void initWRV(FILE * outputFile, double minRelDose, double maxRelDose, long block_size, long block_size_unit_pm){
    fprintf(outputFile, "patdef %ld %ld %ld %f %f 0 0\n", block_size_unit_pm, block_size, block_size, maxRelDose, minRelDose);
    fprintf(outputFile, "vepdef 20 %ld %ld\n", block_size, block_size);
}

void initGTX(FILE *outputFile, double subFieldResolution, long offsetX, long offsetY)
{
     double beamStepSize = 10;
    /**
     * @brief 
        GTX "1.0.0"

        HEADER
        Title                 "Generic Pattern Format, c RAITH nanofabrication"
        Version               "1.44 UPG"
        MainFieldResolution   0.00100000,0.00100000 ! (um)
        SubFieldResolution    0.00050000,0.00050000 ! (um)
        Resolution            10,10 ! 0.01000000,0.01000000 (um) 
        !  BeamStepSize          20,20 ! 0.01000000,0.01000000 (um) 
        BeamStepSize          10,10 ! 0.05, 0.05
        PatternSize           50000,50000 ! 500.00000000,500.00000000 (um) 
        MainFieldSize         50000,50000 ! 500.00000000,500.00000000 (um) 
        SubFieldSize          420,420 ! 4.20000000,4.20000000 (um) 
        HighTension           100000000 ! (mV) 
        SubFieldPlacement     LOWERLEFT
        MainFieldPlacement    MEANDER
        FractureStyle         SUBFIELD HIGHRESOLUTION ! SMOOTHMODE SPIRALMODE
        NrMainFieldBits       20
        NrSubFieldBits        14
        WorkLevel             HIGH
        HeaderOffset          512
        PatternName           "alt_circle_merge.gpf"
        END
     * 
     */
    fprintf(outputFile, "GTX \"1.0.0\"\n\n");
    fprintf(outputFile, "HEADER\n");
    fprintf(outputFile, "Title                 \"Generic Pattern Format, c RAITH nanofabrication\"\n");
    fprintf(outputFile, "Version               \"1.44 UPG\"\n");
    fprintf(outputFile, "MainFieldResolution   0.00100000,0.00100000 ! (um)\n");
    fprintf(outputFile, "SubFieldResolution    %f,%f ! (um)\n", subFieldResolution, subFieldResolution);
    fprintf(outputFile, "Resolution            10,10 ! 0.01000000,0.01000000 (um) \n");
    fprintf(outputFile, "BeamStepSize          %f,%f ! 0.05, 0.05 \n", beamStepSize, beamStepSize);
    fprintf(outputFile, "PatternSize           50000,50000 ! 500.00000000,500.00000000 (um) \n");
    fprintf(outputFile, "MainFieldSize         50000,50000 ! 500.00000000,500.00000000 (um) \n");
    fprintf(outputFile, "SubFieldSize          420,420 ! 4.20000000,4.20000000 (um) \n");
    fprintf(outputFile, "HighTension           100000000 ! (mV) \n");
    fprintf(outputFile, "SubFieldPlacement     LOWERLEFT\n");
    fprintf(outputFile, "MainFieldPlacement    MEANDER\n");
    fprintf(outputFile, "FractureStyle         SUBFIELD HIGHRESOLUTION ! SMOOTHMODE SPIRALMODE\n");
    fprintf(outputFile, "NrMainFieldBits       20\n");
    fprintf(outputFile, "NrSubFieldBits        14\n");
    fprintf(outputFile, "WorkLevel             HIGH\n");
    fprintf(outputFile, "HeaderOffset          512\n");
    fprintf(outputFile, "PatternName           \"zone_plate.gpf\"\n");
    fprintf(outputFile, "END\n\n");

    // For now let's also write the field data
    /**
     * @brief 
FIELD 1,1
  MoveStage 0

  Frequency MSW 0x8000 
  Frequency LSW 0x0

  MSF 10 ! 20
  Main   X LSW 274288 
  Main   Y LSW 274288 
     * 
     */



    fprintf(outputFile, "FIELD 1,1\n");
    fprintf(outputFile, "\tMoveStage 0\n\n");
    fprintf(outputFile, "\tFrequency MSW 0x8000 \n");
    fprintf(outputFile, "\tFrequency LSW 0x0\n\n");
    fprintf(outputFile, "\tMSF %f ! 20\n", beamStepSize);
    fprintf(outputFile, "\tMain   X LSW %ld \n", offsetX);
    fprintf(outputFile, "\tMain   Y LSW %ld \n", offsetY);
}

void initARC(double centerX, double centerY, double nwaUnit, FILE *outputFile)
{
    fprintf(outputFile, "[HEADER]\n");
    fprintf(outputFile, "unit_size %0.2f nm/px\n", nwaUnit * 1000);
    fprintf(outputFile, "Center_um %0.4f %0.4f\n", centerX, centerY);
    fprintf(outputFile, "Center_px %0.2f %0.2f\n", centerX / nwaUnit, centerY / nwaUnit);
    fprintf(outputFile, "\n\n");
    fprintf(outputFile, "[STITCH]\n\n%d, %4d\n%d, %4d\n\n[DATA]\n\n", 0, 0, 0, 0);
}

void initGDS(FILE *outputFile, unsigned char *gdsPost, unsigned char *polyPre, unsigned char *polyPost, unsigned char *polyForm, int layerNumber)
{
    int gdspreamble[102] = {0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                            230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
                            109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
                            243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
                            0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
                            0, 10, 6, 6, 110, 111, 110, 97, 109, 101};

    int gdspostamble[8] = {0, 4, 7, 0, 0, 4, 4, 0};

    // {0, 4, 8, 0, 0, 6, 13, 2, 0, [layer = 1], 0, 6, 14, 2, 0, 0};
    int polypreamble[16] = {0, 4, 8, 0, 0, 6, 13, 2, 0, layerNumber, 0, 6, 14, 2, 0, 0};
    int polypostamble[4] = {0, 4, 17, 0};
    int polyBlockFormat[4] = {0, 44, 16, 3};

    unsigned char gdsPre[102];
    for (int k = 0; k < 102; k++)
        gdsPre[k] = (unsigned char)gdspreamble[k];
    for (int k = 0; k < 8; k++)
        gdsPost[k] = (unsigned char)gdspostamble[k];
    for (int k = 0; k < 16; k++)
        polyPre[k] = (unsigned char)polypreamble[k];
    for (int k = 0; k < 4; k++)
        polyPost[k] = (unsigned char)polypostamble[k];
    for (int k = 0; k < 4; k++)
        polyForm[k] = (unsigned char)polyBlockFormat[k];

    fwrite(gdsPre, sizeof(char), 102, outputFile);
}

void renderGDS(FILE *outputFile, unsigned char *gdsPost)
{
    fwrite(gdsPost, sizeof(char), 8, outputFile);
    fclose(outputFile);
}
void renderARC(FILE *outputFile)
{
    fprintf(outputFile, "[DATA]\n");
    fclose(outputFile);
}
void renderWRV(FILE *outputFile)
{
    fclose(outputFile);
}
void renderGTX(FILE *outputFile)
{
    // fprintf(outputFile, "Activate Merge\n");
    fprintf(outputFile, "END\n");

    fclose(outputFile);
}



// Fractures trapeizoids into two triangles and a horizontal trapezoid for WRV
void fractureAndWriteWRVPoly(double * coords, FILE * outputFile, long clockSpeed, double blockGrid_pixels){
    double heightTol = 1; // shapes with smaller than 0.5 nm height are ignored
    // Need to sort coords
    double x[4] = {coords[0], coords[2], coords[4], coords[6]};
    double y[4] = {coords[1], coords[3], coords[5], coords[7]};
    
    // printf("Presorting: ");
    // printf("[[%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]; [%0.5f,%0.5f]]\n", x[0], y[0],x[1], y[1],x[2], y[2],x[3], y[3]);

    // Sortrows
    sortRows(x, y);
    // printf("Postsorting: ");
    // printf("[[%0.2f,%0.2f]; [%0.2f,%0.2f]; [%0.2f,%0.2f]; [%0.2f,%0.2f]]", x[0], y[0],x[1], y[1],x[2], y[2],x[3], y[3]);
    
    // Find intersections.  Two cases: y[3] -> y[1] (more common), or y[3] -> y[0] (less common)
    double xU, xD;

    // There are two types of trapezoids:
    /**
     *  1) 0 -> 3 is max slope or min slope
     *  2) 0 -> 3 is in between 0->2 and 0->1
     */

    double m03 = atan2((y[3] - y[0]),(x[3] - x[0]));
    double m02 = atan2((y[2] - y[0]),(x[2] - x[0]));
    double m01 = atan2((y[1] - y[0]),(x[1] - x[0]));

    if ( (m03 > m01) == (m02 > m03)){   // rotated rectangle type trap
        // printf(": Used type 1 trapezoid (mid) \n");
        xU = x[3] - (x[3] - x[1]) * (y[3] - y[2]) / (y[3] - y[1]);
        xD = x[2] - (x[2] - x[0]) * (y[2] - y[1]) / (y[2] - y[0]);
    }
    else
    {   // case 2 trapezoid
        //  printf(": Used type 1 trapezoid (L/R) \n");
        // xU = x[3] - (x[3] - x[0])*(y[3] - y[2])/(y[3] - y[0]);
        // xD = x[3] - (x[3] - x[0])*(y[3] - y[1])/(y[3] - y[0]);
        xU = x[0] + (x[3] - x[0]) * (y[2] - y[0]) / (y[3] - y[0]);
        xD = x[0] + (x[3] - x[0]) * (y[1] - y[0]) / (y[3] - y[0]);
    }
    double xDL, xDR, xUL, xUR;
    
    // Find orientation of lower and upper triangles
    if (y[0] == y[1]){ // then set xDL and xDR based on original coords:
        if (x[0] < x[1]){
            xDL = x[0];
            xDR = x[1];
        } else{
            xDL = x[1];
            xDR = x[0];
        }

    } else if (xD < x[1]){
        xDL = xD;
        xDR = x[1];
    }
    else
    {
        xDL = x[1];
        xDR = xD;
    }

    if (y[2] == y[3]){
        if (x[2] < x[3]){
            xUL = x[2];
            xUR = x[3];
        } else{
            xUL = x[3];
            xUR = x[2];
        }
        // printf("Using horizontal top, x[2] = %0.2f, x[3] = %0.2f xUR = %0.2f\n", x[2],x[3], xUR);

    } else if (xU < x[2]){
        xUL = xU;
        xUR = x[2];
    }
    else
    {
        xUL = x[2];
        xUR = xU;
    }

    //     printf("Lower tri: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
    //             (long) x[0], (long)y[0], (long)xDR, (long)y[1], (long)x[0], (long)xDL); // lower tri
    //  printf("Upper tri: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
    //             (long) xUL, (long)y[2], (long)x[3], (long)y[3], (long)xUR, (long)x[3]); // upper tri
    // printf("Middle trap: Trap/%ld %ld %ld %ld %ld %ld %ld\n", clockSpeed,
    //                 (long) xDL, (long)y[1], (long)xUR, (long)y[2], (long)xDR, (long)xUL); // middle trap

    if (y[1] - y[0] > heightTol)
    {
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
    if (y[2] - y[1] > heightTol)
    {
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
               double offsetX, double offsetY, FILE *outputFile)
{
    fprintf(outputFile, "ARC_S %12.4f %10.4f %10.4f %10.4f  %5.2f  %5.2f %5.2f %5.2f D%3.4f S1.00 T-14000\n",
            R / nwaUnit, fmod((theta * 180 / M_PI), 360), dR / nwaUnit, dTheta * 180 / M_PI,
            offsetX / nwaUnit, offsetY / nwaUnit, 0., 0., dose);
}

void exportArcGTX(double R, double dR, double theta, double dTheta, double dose, double block_unit_size_um,
               double offsetX, double offsetY, FILE *outputFile)
{
    
    /**
     * @brief 
        Offset    X LSW 400
        Offset    Y LSW 400

        Circle Diameter Outside LSW 390
        Circle Diameter Inside LSW 100

        Circle Angle Start MSW 0
        Circle Angle Finish MSW 30000
        Activate CIRCLE Pie Ring Merge
     */

    // printf(" offsetX: %f, blokc unit: %f, ratio: %f, cast: %ld", offsetX, block_unit_size_um, offsetX/block_unit_size_um, (long) (offsetX/block_unit_size_um));

    fprintf(outputFile, "\tOffset    X LSW %ld\n", (long)(offsetX / block_unit_size_um));
    fprintf(outputFile, "\tOffset    Y LSW %ld\n", (long)(offsetY / block_unit_size_um));
    fprintf(outputFile, "\tCircle Diameter Outside LSW %ld\n", (long)(2*(R + dR) / block_unit_size_um));
    fprintf(outputFile, "\tCircle Diameter Inside LSW %ld\n", (long)(2*R / block_unit_size_um));
    fprintf(outputFile, "\tCircle Angle Start MSW %ld\n", (long)(1000 * fmod(theta * 180 / M_PI, 360))); // milidegrees
    fprintf(outputFile, "\tCircle Angle Finish MSW %ld\n", (long)(1000 * fmod(180 / M_PI * (theta + dTheta), 360))); // milidegrees
    // fprintf(outputFile, "Activate CIRCLE Pie Ring Merge\n\n");
    fprintf(outputFile, "\tActivate CIRCLE Pie Ring\n\n");


}

void encode32(long aCoord, char *cPart)
{
    cPart[0] = (aCoord >> 24) & 255;
    cPart[1] = (aCoord >> 16) & 255;
    cPart[2] = (aCoord >> 8) & 255;
    cPart[3] = (aCoord)&255;
}

void encodePoly32(double * coords, char * cCoords){
    char cPart[4];
    for (int k = 0; k < 10; k++){
        encode32((long) coords[k], cPart);
        for (int m = 0; m < 4; m++){
            cCoords[k*4 + m] = cPart[m];
        }
    }
}



void exportPolygon(double * coords, unsigned char * polyPre, unsigned char * polyPost,
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
