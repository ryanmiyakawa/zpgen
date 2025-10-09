#ifndef PATTERN_FILE_UTILS_H
#define PATTERN_FILE_UTILS_H

#include <stdio.h>

// Function Prototypes

double getDoseFromBias(double dr, double bias);

void writeSupportTxtHeader(double dR1, double dRN, double bias, FILE *supportTextFile);

void writeSupportTxtZoneDose(long clockSpeed, int n, double Rn, FILE *supportTextFile);

void initWRV(FILE * outputFile, double minRelDose, double maxRelDose, long block_size, long block_size_unit_pm);

void initGTX(FILE *outputFile, double subFieldResolution, long offsetX, long offsetY);

void initARC(double centerX, double centerY, double nwaUnit, FILE *outputFile);

void initGDS(FILE *outputFile, unsigned char *gdsPost, unsigned char *polyPre, unsigned char *polyPost, unsigned char *polyForm, int layerNumber);

void renderGDS(FILE *outputFile, unsigned char *gdsPost);

void renderARC(FILE *outputFile);

void renderWRV(FILE *outputFile);

void renderGTX(FILE *outputFile);

void fractureAndWriteWRVPoly(double * coords, FILE * outputFile, long clockSpeed, double blockGrid_pixels);

void exportArc(double R, double dR, double theta, double dTheta, double dose, double nwaUnit,
               double offsetX, double offsetY, FILE *outputFile);
void exportArcGTX(double R, double dR, double theta, double dTheta, double dose, double block_unit_size_um,
               double offsetX, double offsetY, FILE *outputFile);

void encode32(long aCoord, char *cPart);

void encodePoly32(double * coords, char * cCoords, int numCoords);


void exportPolygon(double * coords, unsigned char * polyPre, unsigned char * polyPost,
                   unsigned char * polyForm, FILE * outputFile, int File_format, long clockSpeed, double blockGrid_pixels, int numVertices);


#endif // PATTERN_FILE_UTILS_H