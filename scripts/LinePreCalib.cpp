#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
//#include <TTree.h>
#include "TH1F.h"
//#include <TH2F.h>
#include <TH1D.h>
//#include <TH2D.h>
#include "TF1.h"
#include "TLegend.h"
//#include <TCut.h>
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
//#include <TGraph.h>
//#include <TGraphErrors.h>
//#include <TIterator.h>
//#include <TList.h>
//#include <TMultiGraph.h>
#include "TMath.h"
//#include "TApplication.h"
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>
//#include <TLatex.h>
//#include <TTimeStamp.h>
//#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
//#include <TFitResult.h>

//#include "../source/include/LineStruct.h"
#include "../source/include/GatorCalibFunc.h"
#include "../source/include/GatorLoadData.h"
#include "../source/include/GatorCalibFitters.h"


using namespace std;


void LinePreCalib(){
	
	const string archivedir("/Users/francesco/PhD/Gator/data_analysis/calibrations/archive/");
	string calibset = string("23-12-2013/54Mn");
	
	string datadir = archivedir + calibset + string("/");
	
	bool dofit  = false;
	
	double MCAlowCh = 1400;
	double MCAupCh = 1420;
	double initMean = 7873;
	double initAmpl = 1.5e2;
	double initTail = 1.0;
	double initSigma = 5.0;
	double initBeta = 1e-1;
	double initStep = 1.0;
	double initCost = 0.1;
	
	
	//Load the spectrum of the calibs
	double calibtime = 0.;
	TH1F* MCAhisto = loadSpe(datadir.c_str(),calibtime);
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	//Make a copy of the ROI of the original histo and put it in a new histogram
	int firstbin = MCAhisto->FindBin(MCAlowCh);
	int lastbin = MCAhisto->FindBin(MCAupCh);
	double xmin = MCAhisto->GetBinLowEdge(firstbin);
	double xmax = MCAhisto->GetBinLowEdge(lastbin+1);
	int nbins = 1 + lastbin - firstbin;
	
	
	TH1D* tmphisto = new TH1D("tmphisto",";MCA channel;Counts",nbins,xmin,xmax);
	
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = MCAhisto->GetBinCenter(bin);
		tmphisto->Fill(bincenter,MCAhisto->GetBinContent(bin));
	}
	
	TCanvas *c1 = new TCanvas("c1");
	
	
	
	
	//Fitting part with MINUIT
	
	
	//Start with the fitting process
	if(dofit){
		
		tmphisto->Draw();
		
		TF1* ff_MCA = new TF1("ff_MCA",peakFitFunc,MCAlowCh,MCAupCh,7);
		ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","StepAmpl","Constant");
		ff_MCA->SetLineColor(kRed);
		ff_MCA->SetLineWidth(2);
	
		ff_MCA -> SetParameter("Mean",initMean);
		ff_MCA -> SetParameter("Ampl",initAmpl);
		ff_MCA -> SetParameter("Tail",initTail);
		ff_MCA -> SetParameter("Sigma",initSigma);
		ff_MCA -> SetParameter("Beta",initBeta);
		ff_MCA -> SetParameter("StepAmpl",initStep);
		ff_MCA -> SetParameter("Constant",initCost);
	
	
		//ff_MCA -> SetParLimits(0,line.MCAlowch,line.MCAupch);
		//ff_MCA -> SetParLimits(1,0.,3*line.ampl);
		//ff_MCA -> SetParLimits(2,0.,1.);
		//ff_MCA -> SetParLimits(3,line.sigma/4.,4*line.sigma);
		ff_MCA -> SetParLimits(4,0.,1e3);
		//ff_MCA -> SetParLimits(5,0.,5*line.step);
		//ff_MCA -> SetParLimits(1,0.,2*line.cost);
	
		tmphisto->Fit(ff_MCA,"LM");
	}else{
		MCAhisto->Draw();
		tmphisto->SetLineColor(kRed);
		tmphisto->Draw("same");
	}
	
	return;
}

#include "../source/src/loadSpe.cc"
#include "../source/src/GatorCalibFitters.cc"