#ifndef GATORCALIBFITTERS_H
#define GATORCALIBFITTERS_H

#include "TH1D.h"
#include "GatorStructs.h"

Double_t peakFitFuncA(Double_t* x,Double_t* par);
Double_t peakFitFuncB(Double_t* x,Double_t* par);

void amplInit(TH1D* histo, CalibLine& line);
void costInit(TH1D* histo, CalibLine& line);
void stepInit(TH1D* histo, CalibLine& line);
void sigmaInit(TH1D* histo, CalibLine& line);

Double_t resolFunc (Double_t* x,Double_t* par);


#endif
