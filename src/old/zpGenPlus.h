
#ifndef ZPGENPLUS_H
#define ZPGENPLUS_H

void makeZP(double zTol, double lambda_nm, double p, double q, double obscurationSigma, 
    double NA, int nZerns, double *orders, int customMaskIdx, double beta, double CRAAz, 
    double CRA, double anamorphicFac, double ZPCPhase, double APD, double APD_window, 
    double ZPCR2, double ZPCR1, double bias_nm, int File_format, int Opposite_Tone, 
    int FSIdx, double buttressGapWidth, double buttressPeriod, int setToCenter, 
    long block_size, int NoP, int IoP, double blockGrid_pm, int layerNumber, int nwaUnitSelection,
    char * fileName);


#endif
