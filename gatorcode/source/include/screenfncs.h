#ifndef SCREENFNCS_H
#define SCREENFNCS_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
//#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>
//#include <TCut.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
//#include <TIterator.h>
//#include <TList.h>
//#include <TMultiGraph.h>
#include <TMath.h>
#include <TApplication.h>
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>
//#include <TLatex.h>
//#include <TTimeStamp.h>
#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
#include <TFitResult.h>


using namespace std;

//Only for this functin the definition is made in the header
bool fexist(string filename);

TH1F* loadSPE(const char* dir, double& aqtime);
TH1F* loadSingleSPE(const char* name, double& aqtime);
TH1F* loadSpe(const char* dir, double& aqtime);

string formatdigits1(double var, double var_err, int dig = 2);

string formatdigits2(double var, int dig = 3);

TH1F* convert_histo_ENR(TH1F* MCAspect,const char* calibdir);


#endif