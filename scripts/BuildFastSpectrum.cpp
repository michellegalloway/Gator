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

TH1D* loadSPE(const char* dir, Double_t& aqtime);
TH1D* convert_histo_ENR(TH1D* hADC, const char* calibdir);
TCanvas* c1;

void BuildFastSpectrum(){

	//Style

	gStyle  ->SetOptStat(0);
    gStyle  ->SetOptFit(0);

    //------ define color gradient
    const Int_t NRGBs = 5;
    const Int_t NCont = 10;

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
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadTopMargin(0.01);
	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadLeftMargin(0.05);
	//Finish of ATP common makeup
	
	
	string workdir = "/Users/francesco/PhD/Gator/analysis/samples/TiGr1_Uniti/";
	string outfilename = workdir + string("TiGr1_Uniti_spect.root");
	
	
	//string SPEdir = "/Users/francesco/PhD/Gator/analysis/sessions/XeActivation/SPE/activ/";
	//string SPEdir = "/Users/francesco/PhD/Gator/analysis/samples/MgF2_PMT/SPE/";
	//string calibdir = "/Users/francesco/PhD/Gator/calibrations/archive/2015.08.07/";
	
	//string SPEdir = "/Users/francesco/PhD/Gator/analysis/samples/PTFE_40kg/SPE/";
	string SPEdir = "/Users/francesco/PhD/Gator/analysis/samples/PTFE_40kg/SPEreduced/";
	string calibdir = "/Users/francesco/PhD/Gator/calibrations/archive/2013.12.23/";
	
	//string SPEdir = "/Users/francesco/PhD/Gator/analysis/samples/TiGr1_Uniti/SPE/";
	//string calibdir = "/Users/francesco/PhD/Gator/calibrations/archive/2012.06.19/";
	
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/analysis/sessions/XeActivation/SPE/pre_activ/";
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/analysis/sessions/BG_137Cs/EmptyCavity/";
	//string calibdir_bg = "/Users/francesco/PhD/Gator/calibrations/archive/2013.12.23/";
	
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/background/archive/2011/";
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/background/archive/2012/";
	//string calibdir_bg = "/Users/francesco/PhD/Gator/calibrations/archive/2012.06.19/";
	
	string backgroundSPEdir = "/Users/francesco/PhD/Gator/background/archive/2014/";
	string calibdir_bg = "/Users/francesco/PhD/Gator/calibrations/archive/2013.12.23/";
	
	
	int rebin = 8;
	
	double minEn = 5;
	double maxEn = 2700;
	
	string sampleLeg = "Sample";
	string bgLeg = "Background";
	
	//string sampleLeg = "Xenon bottle - post activation";
	//string bgLeg = "Xenon bottle - pre activation";
	
	Double_t aqtime, aqtime_bg;
	
	TH1D* hADC = loadSPE(SPEdir.c_str(), aqtime);
	hADC->SetName("hSampleMCA");
	
	TH1D* hADC_bg = loadSPE(backgroundSPEdir.c_str(), aqtime_bg);
	hADC_bg->SetName("hBgMCA");
	
	
	c1 = new TCanvas("c1","",1200,400);
	c1->SetLogy();
	
	TH1D* histoENR = convert_histo_ENR(hADC, calibdir.c_str());
	histoENR -> SetName("histoENR");
	histoENR -> SetTitle("");
	if(rebin>1) histoENR -> Rebin(rebin);
	histoENR -> Scale(1./(aqtime/(24*3600)),"width");
	histoENR -> SetLineWidth(2);
	histoENR -> SetLineColor(kRed-3);
	//histoENR -> SetLineColor(kRed);
	histoENR -> GetXaxis() -> SetRangeUser(minEn,maxEn);
	histoENR -> GetXaxis() -> SetTitle("Energy [keV]");
	histoENR -> GetYaxis() -> SetTitle("Differential rate [counts keV^{-1} day^{-1}]");
	
	//MAKEUP
	histoENR->GetXaxis()->SetLabelFont(132); //Times New Roman
	histoENR->GetXaxis()->SetLabelSize(0.045);
	histoENR->GetXaxis()->SetTitleFont(132); //Times New Roman
	histoENR->GetXaxis()->SetTitleSize(0.045);
	histoENR->GetXaxis()->SetTitleOffset(1.1);
	
	histoENR->GetYaxis()->SetLabelFont(132); //Times New Roman
	histoENR->GetYaxis()->SetLabelSize(0.045);
	histoENR->GetYaxis()->SetTitleFont(132); //Times New Roman
	histoENR->GetYaxis()->SetTitleSize(0.045);
	histoENR->GetYaxis()->SetTitleOffset(0.55);
	
	
	
	
	TH1D* bg_histoENR = convert_histo_ENR(hADC_bg,calibdir_bg.c_str());
	bg_histoENR -> SetName("bg_histoENR");
	bg_histoENR -> SetTitle("");
	if(rebin>1) bg_histoENR -> Rebin(rebin);
	bg_histoENR -> Scale(1./(aqtime_bg/(24*3600)),"width");
	bg_histoENR -> SetLineWidth(1);
	bg_histoENR -> SetLineColor(kGray+2);
	//bg_histoENR -> SetLineColor(kBlue);
	bg_histoENR -> GetXaxis() -> SetRangeUser(minEn,maxEn);
	bg_histoENR -> GetXaxis() -> SetTitle("Energy [keV]");
	bg_histoENR -> GetYaxis() -> SetTitle("Differential rate [counts keV^{-1} day^{-1}]");
	
	//MAKEUP
	bg_histoENR->GetXaxis()->SetLabelFont(132); //Times New Roman
	bg_histoENR->GetXaxis()->SetLabelSize(0.045);
	bg_histoENR->GetXaxis()->SetTitleFont(132); //Times New Roman
	bg_histoENR->GetXaxis()->SetTitleSize(0.045);
	bg_histoENR->GetXaxis()->SetTitleOffset(1.1);
	
	bg_histoENR->GetYaxis()->SetLabelFont(132); //Times New Roman
	bg_histoENR->GetYaxis()->SetLabelSize(0.045);
	bg_histoENR->GetYaxis()->SetTitleFont(132); //Times New Roman
	bg_histoENR->GetYaxis()->SetTitleSize(0.045);
	bg_histoENR->GetYaxis()->SetTitleOffset(0.55);
	
	histoENR -> Draw();
	bg_histoENR -> Draw("same");
	
	
	//return;
	
	TLegend* leg = new TLegend(0.60,0.90,0.90,0.80);
	leg -> AddEntry(histoENR,sampleLeg.c_str(),"L");
	leg -> AddEntry(bg_histoENR,bgLeg.c_str(),"L");
	leg -> SetFillColor(kWhite);
	leg -> SetBorderSize(0);
	leg -> SetTextFont(62);
	leg -> Draw();
	
	cout << "\nAcquisition time sample: " << aqtime << " s ---> " << aqtime/(24*3600) << " d" << endl;
	cout << "Acquisition time bckg: " << aqtime_bg << " s ---> " << aqtime_bg/(24*3600) << " d" << endl;
	
	//c1->SaveAs(outfilename.c_str());
	
	return;		
}

#include "../source/src/loadSPE.cc"
#include "../source/src/convert_histo_ENR.cc"