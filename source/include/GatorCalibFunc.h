#ifndef GATORCALIB_H
#define GATORCALIB_H

#include <vector>
#include <string>

#include "GatorStructs.h"

#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

vector<string>& loadConfFile(string filename);

CalibLine lineInit(string line);

bool doFitBAT(TH1D* MCAhisto, CalibLine& line);

bool doFitLL(TH1D* MCAhisto, CalibLine& line);

vector<CalibLine> LoadLinesFromTree(string datadir);

TGraphErrors* plotResiduals(TF1* calibfunc,TGraphErrors* grPlot,TFitResultPtr& resFit);

#endif
