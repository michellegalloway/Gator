#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include <TStyle.h>
#include "TFile.h"
//#include <TTree.h>
#include "TH1F.h"
//#include <TH2F.h>
#include <TH1D.h>
//#include <TH2D.h>
#include "TParameter.h"
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
#include "TApplication.h"
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


using namespace std;

TH1F* loadSPE(const char* dir, Double_t& aqtime);
TH1F* convert_histo_ENR(TH1F* hADC, const char* calibdir);

void BgFastSpectrum(){

	//Style

	gStyle  ->SetOptStat(0);
    gStyle  ->SetOptFit(0);

    //------ define color gradient
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    //gStyle->SetNumberContours(30);
    //--------------------------------------------------------------------------
	gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(10);
	gStyle->SetStatColor(10);
	gStyle->SetStatFont(42);
	//Finish of ATP common makeup
	
	
	const Double_t hour = 3600; //Seconds in an hour
	const Double_t day = 24*hour; //Seconds in a day
	const Double_t Mass = 2.2; //Mass of the sensitive Gator cristal in kg
	
	
	string bgdir = "/Users/francesco/PhD/Gator/analysis/background/archive/2015/";
	string calibdir = "/Users/francesco/PhD/Gator/analysis/calibrations/archive/2015.08.07/";
	
	string outfilename("/Users/francesco/Phd/Gator/analysis/background/archive/2015/bg2015_spec.root");
	
	string sampleLeg("");
	
	Double_t aqtime;
	
	TH1F* hMCA = loadSPE(bgdir.c_str(), aqtime);
	hMCA->SetNameTitle("bgMCA","; MCA channel; Counts");
	
	TH1F* histoENR = convert_histo_ENR(hMCA, calibdir.c_str());
	
	cout << "\nAcquisition time: " << aqtime/day << " days" << endl;
	
	TParameter<double>* tp_aqtime = new TParameter<double>("aqtime",aqtime);
	
	//hMCA -> Scale(1./(aqtime/day));
	
	histoENR -> SetName("bgENR2");
	histoENR -> SetTitle("; Energy [keV]; Counts keV^{-1} day^{-1}");
	histoENR -> Rebin(8);
	histoENR -> Scale(1./(aqtime/day),"width");
	histoENR -> SetLineWidth(2);
	histoENR -> SetLineColor(kRed);
	
	
	TCanvas* c1 = new TCanvas("c1","Background 2015 (MCA)");
	c1->SetLogy();
	hMCA -> Draw();
	
	TCanvas* c2 = new TCanvas("c2","Background 2015 (ENR)");
	c2->SetLogy();
	histoENR -> Draw();
	
	TFile* outFile = new TFile((outfilename).c_str(),"update");
	outFile->WriteTObject(c1);
	outFile->WriteTObject(c2);
	outFile->WriteTObject(hMCA);
	outFile->WriteTObject(histoENR,0,"overwrite");
	outFile->WriteTObject(tp_aqtime);
	//outFile->Close();
	
	return;
}

//Here the libraries
#include "/Users/francesco/PhD/Gator/source/src/loadSPE.cc"
#include "/Users/francesco/PhD/Gator/source/src/convert_histo_ENR.cc"