#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include <TTree.h>
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
#include <TGraph.h>
#include <TGraphErrors.h>
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
//#include "../source/include/GatorCalibFunc.h"
//#include "../source/include/GatorLoadData.h"
//#include "../source/include/GatorCalibFitters.h"
#include "../source/include/screenfncs.h"



using namespace std;


void ShowCalib(){
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1100);
	gStyle->SetStatBorderSize(0);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(10);
	gStyle->SetStatColor(10);
	gStyle->SetStatFont(42);
	gStyle->SetMarkerStyle(7);
	gStyle->SetMarkerColor(2);
	
	const string archivedir("/Users/francesco/PhD/Gator/calibrations/archive/");
	string calibset = string("2015.07.31/137Cs");
	
	string datadir = archivedir + calibset + string("/");
	
	//Load the spectrum of the calibs
	double calibtime;
	TH1F* MCAhisto = loadSpe(datadir.c_str(),calibtime);
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	TCanvas *c1 = new TCanvas("c1");
	
	MCAhisto -> Draw();
	
	return;
	
	//CALIBRATION PLOT
	TCanvas *c2 = new TCanvas("c2","c2");
	
	TFile *linesfile = TFile::Open( (datadir+string("fittedlines.root")).c_str(), "read" );
	
	TTree *linestree = (TTree*)linesfile->Get("linestree");
	
	int nLines = linestree->GetEntries();
	
	double litEn, litEn_err, mean, mean_err, sigma, sigma_err;
	
	linestree->SetBranchAddress("litEn", &litEn);
	linestree->SetBranchAddress("litEn_err", &litEn_err);
	linestree->SetBranchAddress("mean", &mean);
	linestree->SetBranchAddress("mean_err", &mean_err);
	linestree->SetBranchAddress("sigma", &sigma);
	linestree->SetBranchAddress("sigma_err", &sigma_err);
	
	TGraphErrors *calibplot = new TGraphErrors();
	for(int iEv=0; iEv<nLines; iEv++){
		
		linestree->GetEntry(iEv);
		
		calibplot->SetPoint(iEv, mean, litEn);
		calibplot->SetPointError(iEv, mean_err, litEn_err);
		
	}
	//calibplot->SetMarkerStyle(20);
	calibplot->Draw("AP");
	
	
	TF1 *calib_fcn_fw = new TF1("calib_fcn_fw", "pol2", 0, 16384);
	calib_fcn_fw->SetLineWidth(2);
	calibplot->Fit(calib_fcn_fw);
	//calib_fcn_fw->Draw("SAME");
	
	
	//RESIDUAL PLOT
	TCanvas *c3 = new TCanvas("c3","c3");
	
	TFile *calibration = TFile::Open( (datadir+string("calibration.root")).c_str(), "read" );
	
	TGraphErrors *gr_cal_res_fw = (TGraphErrors*)calibration->Get("gr_cal_res_fw");
	gr_cal_res_fw->SetMarkerStyle(20);
	gr_cal_res_fw->SetMarkerColor(kBlack);
	gr_cal_res_fw->SetLineColor(kBlack);
	gr_cal_res_fw->SetTitle("");
	gr_cal_res_fw->Draw("AP");
	
	TF1 *zeroline = new TF1("zeroline", "pol0", 0, 2800);
	zeroline->SetParameter(0,0);
	zeroline->SetLineWidth(2);
	zeroline->SetLineColor(kRed);
	
	zeroline->Draw("SAME");
	
	
	//ENERGY RESOLUTION PLOT
	TCanvas *c4 = new TCanvas("c4","c4");
	
	double p0 = calib_fcn_fw->GetParameter(0);
	double p0err = calib_fcn_fw->GetParError(0);
	double p1 = calib_fcn_fw->GetParameter(1);
	double p1err = calib_fcn_fw->GetParError(1);
	double p2 = calib_fcn_fw->GetParameter(2);
	double p2err = calib_fcn_fw->GetParError(2);
	
	double sigmaEn, sigmaEn_err;
	
	TGraphErrors *sigmaENRplot = new TGraphErrors();
	for(int iEv=0; iEv<nLines; iEv++){
		
		linestree->GetEntry(iEv);
		
		sigmaEn = (p1+4*p2*mean)*sigma;
		sigmaEn_err = sqrt( pow(p1err*sigmaEn,2) + pow(4*mean*sigmaEn*p2err,2) + pow(4*p2*sigma*mean_err,2) + pow((p1+4*p2*mean)*mean_err,2) );
		
		sigmaENRplot->SetPoint(iEv, litEn, sigmaEn);
		sigmaENRplot->SetPointError(iEv, litEn_err, sigmaEn_err);
		
	}
	sigmaENRplot->SetTitle("; Energy [keV]; #sigma_{E} [keV]");
	sigmaENRplot->GetXaxis()->CenterTitle();
	sigmaENRplot->SetMarkerColor(kBlack);
	sigmaENRplot->SetMarkerStyle(20);
	sigmaENRplot->Draw("AP");
	
	TF1 *resol_enr_fnc = new TF1("resol_enr_fnc","sqrt([0]+[1]*x+[2]*x*x)",50,2750);
	
	resol_enr_fnc->SetParameters(1e-1,5e-4,0);
	resol_enr_fnc->SetLineColor(kRed);
	resol_enr_fnc->SetLineWidth(2);
	sigmaENRplot->Fit(resol_enr_fnc);
	
	return;
	
}

#include "../source/src/loadSPE.cc"
//#include "../source/src/GatorCalibFitters.cc"