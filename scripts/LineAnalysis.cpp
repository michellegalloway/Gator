#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
//#include <TStyle.h>
#include "TFile.h"
//#include <TTree.h>
#include "TH1F.h"
//#include <TH2F.h>
//#include <TH1D.h>
//#include <TH2D.h>
#include "TF1.h"
#include "TLegend.h"
//#include <TCut.h>
#include "TAxis.h"
#include "TCanvas.h"
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
#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
//#include <TFitResult.h>


using namespace std;

//Declare global variables (can be used in the root interactive session after the script has run)
Double_t lineADC, sigmaADC;
TH1F* hADC;


TH1F* loadSPE(const char* dir, Double_t& aqtime);
//void linecounts2(const TH1F*, TCanvas*, const Double_t, const Double_t, const Double_t, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&, Bool_t&, Double_t&, const TH1F* bg_histoADC = 0, Double_t bg_aqtime = 0, Double_t WidthScale = 3.0);


//Function definitions inclusion
#include "/Users/francesco/PhD/Gator/data_analysis/samples/src/loadSPE.cpp"



void LineAnalysis(){

	//HARDCODED STUFF
	string workdir = "/Users/francesco/PhD/Gator/data_analysis/samples/PMTs/R11410-21/CeramicStems/";
	
	string SPEdir = "/Users/francesco/PhD/Gator/data_analysis/samples/PMTs/R11410-21/CeramicStems/SPE/";
	
	string calibdir = "/Users/francesco/PhD/Gator/data_analysis/calibrations/archive/04-10-2012/";
	
	string backgroundSPEdir = "/Users/francesco/PhD/Gator/data_analysis/background/";
	
	string calibdir_bg = "/Users/francesco/PhD/Gator/data_analysis/calibrations/archive/19-06-2012/";

	string outfilename = workdir + string("NironitSS_spectrum.C");
	
	Double_t lineENR = 1764.49;
	Double_t BRxEff = 0.00146;
	//FINISHED HARDCODED STUFF
	
	string calfilename = calibdir + string("calibration.root");
	string calfilename_bg = calibdir_bg + string("calibration.root");
	
	
	//Load sample and background raw spectra
	Double_t aqtime, aqtime_bg;
	hADC = loadSPE(SPEdir.c_str(), aqtime);
	TH1F* hADC_bg = loadSPE(backgroundSPEdir.c_str(), aqtime_bg);
	
	
	//Load calibrations
	TFile* calfile = new TFile(calfilename.c_str(),"read");
	TF1* calib_fcn_bw = (TF1*)calfile -> Get("calib_fcn_bw");
	TF1* resol_ADC_func = (TF1*)calfile -> Get("resol_ADC_func");
	calfile -> Close();
	
	
	lineADC = calib_fcn_bw->Eval(lineENR);
	sigmaADC = resol_ADC_func->Eval(lineADC);
	
	cout << endl << "lineADC: " << lineADC << endl;
	cout << "sigmaADC: " << sigmaADC << endl;
	
	
	//Draw the spectra and the counting regions
	TCanvas* c1 = new TCanvas("c1","");
	c1->SetLogy();
	
	hADC->SetLineColor(kRed);
	hADC->Draw();
	
	
	Double_t y_linelow = hADC->GetBinContent(hADC->FindBin(lineADC))*TMath::Exp(-0.5);
	Double_t y_lineup = hADC->GetBinContent(hADC->FindBin(lineADC))*TMath::Exp(0.5);
	
	TLine* lineC = new TLine(lineADC,y_linelow,lineADC,y_lineup);
	lineC->SetLineColor(kBlack);
	lineC->SetLineWidth(2);
	lineC->Draw("same");
	
	TLine* lineR = new TLine(lineADC+3.*sigmaADC,y_linelow,lineADC+3.*sigmaADC,y_lineup);
	lineR->SetLineColor(kRed);
	lineR->SetLineWidth(2);
	lineR->Draw();
	
	TLine* lineL = new TLine(lineADC-3.*sigmaADC,y_linelow,lineADC-3.*sigmaADC,y_lineup);
	lineL->SetLineColor(kRed);
	lineL->SetLineWidth(2);
	lineL->Draw();
	
	TLine* lineRR = new TLine(lineADC+6.*sigmaADC,y_linelow,lineADC+6.*sigmaADC,y_lineup);
	lineRR->SetLineColor(kBlue);
	lineRR->SetLineWidth(2);
	lineRR->Draw();
	
	TLine* lineLL = new TLine(lineADC-6.*sigmaADC,y_linelow,lineADC-6.*sigmaADC,y_lineup);
	lineLL->SetLineColor(kBlue);
	lineLL->SetLineWidth(2);
	lineLL->Draw();
	
	c1->Modified();
	
	return;
}

