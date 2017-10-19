#include "GatorGlobals.hh"
#include "GatorStructs.h"
#include "GatorCalibFunc.h"
//#include "GatorCalibFitters.h"
//#include "GatorCalibScr.h"
//#include "GatorDataLoader.hh"

#include <string>
#include <cstdlib>
#include <unistd.h>

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TFile.h"

#include <BAT/BCAux.h>
//#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>


using namespace std;

void usage();

TApplication* theApp;

//Load the data from the root file of the fitted lines
//const string archivedir("/home/atp/fpiastra/Gator/calibrations/archive/");
const string archivedir("/Users/francesco/PhD/Gator/calibrations/archive/");

int main(int argc, char* argv[]){
	
	theApp = new TApplication("App",0,0);
	
	TCanvas *c1 = new TCanvas("c1","MCA -> Energy conversion");
	c1 -> Divide(1,2);
	
	TCanvas *c2 = new TCanvas("c2","Energy -> MCA conversion");
	c2 -> Divide(1,2);
	
	TCanvas *c3 = new TCanvas("c3","MCA resolution");
	
	
/*	
//The following lines are for closing the application when the canvas c1 is closed.
#ifndef __CINT__
	c1 -> Connect("TCanvas", "Closed()", "TApplication", gApplication,
				  "Terminate()");
#endif
	
/*
	int mickey = 0;
	char *mouse[10];
	theApp = new TApplication("theApp", &mickey, mouse );
*/
	
	char c;
	
	bool calibset_flag = false;
	
	string calibset("");
	
    // parse switches
    while((c = getopt(argc,argv,"d:h")) != -1)
    {
      switch(c)	{
        
        case 'd': 
          calibset_flag = true;
          calibset = optarg;
          break;
        
        case 'h':
			usage();
			return 0;
          	break;
        
        default:
		cout << endl << argv[0] << ":ERROR --> \"-" << c << " doesn't match with any option!!!" << endl << endl;
          usage();
		  return 0;
      }
    }
	
	if(!calibset_flag){
		cout << endl << argv[0] << ":ERROR --> The calibration set is missing. Exiting!!!" << endl << endl;
		usage();
		return(-1);
	}
	
	
	
	
	vector<CalibLine> linesVec = LoadLinesFromTree(archivedir+calibset+string("/"));
	
	if(linesVec.size()==0){
		cout << "ERROR ::::: No entries loaded in the vector of lines.\nExiting!!!" << endl;
		exit(-1);
	}
	
	int nLines = linesVec.size();
	
	//Plotting and fitting the calibration and residuals
	
	TGraphErrors* gr_cal_fw = new TGraphErrors(nLines);
	gr_cal_fw -> GetXaxis()->SetTitle("MCA channel");
	gr_cal_fw -> GetYaxis()->SetTitle("Energy [KeV]");
	gr_cal_fw -> SetMarkerStyle(7);
	gr_cal_fw -> SetMarkerSize(2);
	
	TGraphErrors* gr_cal_bw = new TGraphErrors(nLines);
	gr_cal_bw -> GetXaxis()->SetTitle("Energy [KeV]");
	gr_cal_bw -> GetYaxis()->SetTitle("MCA channel");
	gr_cal_bw -> SetMarkerStyle(7);
	gr_cal_bw -> SetMarkerSize(2);
	
	TGraphErrors* gr_resol_ADC = new TGraphErrors(nLines);
	gr_resol_ADC -> GetXaxis()->SetTitle("MCA channel");
	gr_resol_ADC -> GetYaxis()->SetTitle("#sigma MCA");
	gr_resol_ADC -> SetMarkerStyle(7);
	gr_resol_ADC -> SetMarkerSize(2);
	
	TGraphErrors* gr_cal_res_fw = NULL; //Fill this with the the function plotResiduals
	TGraphErrors* gr_cal_res_bw = NULL; //Fill this with the the function plotResiduals
	
	
	int point =0;
	
	for(int iLine=0; iLine<nLines; iLine++){
		gr_cal_fw -> SetPoint(iLine,linesVec.at(iLine).mean,linesVec.at(iLine).litEn);
		gr_cal_fw -> SetPointError(iLine,linesVec.at(iLine).mean_err,linesVec.at(iLine).litEn_err);
		
		gr_cal_bw -> SetPoint(iLine,linesVec.at(iLine).litEn,linesVec.at(iLine).mean);
		gr_cal_bw -> SetPointError(iLine,linesVec.at(iLine).litEn_err,linesVec.at(iLine).mean_err);
		
		gr_resol_ADC -> SetPoint(iLine,linesVec.at(iLine).mean,linesVec.at(iLine).sigma);
		gr_resol_ADC -> SetPointError(iLine,linesVec.at(iLine).mean_err,linesVec.at(iLine).sigma_err);
		
	}
	
	
	
	//Draw the plots and do the fits
	c1 -> cd(1);
	gr_cal_fw -> Draw("AP");
	
	TF1* calib_fcn_fw = new TF1("calib_fcn_fw","pol2");
	calib_fcn_fw -> SetRange(0., 20000);
	calib_fcn_fw -> SetLineColor(2);
	calib_fcn_fw -> SetLineWidth(1);
	TFitResultPtr resFit_fw = gr_cal_fw -> Fit(calib_fcn_fw,"S");
	
	
	c2 -> cd(1);
	gr_cal_bw -> Draw("AP");
	
	TF1* calib_fcn_bw = new TF1("calib_fcn_bw","pol2");
	calib_fcn_fw -> SetRange(0., 3000);
	calib_fcn_bw -> SetLineColor(2);
	calib_fcn_bw -> SetLineWidth(1);
	TFitResultPtr resFit_bw = gr_cal_bw -> Fit(calib_fcn_bw,"S");
	
	
	//Build the plots for the residuals
	gr_cal_res_fw = plotResiduals(calib_fcn_fw,gr_cal_fw,resFit_fw);
	gr_cal_res_fw -> SetMarkerStyle(7);
	gr_cal_res_fw -> SetMarkerSize(2);
	
	gr_cal_res_bw = plotResiduals(calib_fcn_bw,gr_cal_bw,resFit_bw);
	gr_cal_res_bw -> SetMarkerStyle(7);
	gr_cal_res_bw -> SetMarkerSize(2);
	
	
	c1 -> cd(2);
	gr_cal_res_fw -> GetXaxis() -> SetTitle("MCA channel");
	gr_cal_res_fw -> GetYaxis() -> SetTitle("Energy reconstruction residuals [keV]");
	gr_cal_res_fw -> Draw("AP");
	
	TLine* zeroline_fw = new TLine(gr_cal_res_fw->GetXaxis()->GetXmin(),0,gr_cal_res_fw->GetXaxis()->GetXmax(),0);
	zeroline_fw -> SetLineColor(2);
	zeroline_fw -> SetLineWidth(2);
	zeroline_fw -> Draw("same");
	
	c1 -> Update();
	
	
	c2 -> cd(2);
	gr_cal_res_bw -> GetXaxis() -> SetTitle("Energy [keV]");
	gr_cal_res_bw -> GetYaxis() -> SetTitle("MCA channels reconstruction residuals");
	gr_cal_res_bw -> Draw("AP");
	
	TLine* zeroline_bw = new TLine(gr_cal_res_bw->GetXaxis()->GetXmin(),0,gr_cal_res_bw->GetXaxis()->GetXmax(),0);
	zeroline_bw -> SetLineColor(2);
	zeroline_bw -> SetLineWidth(2);
	zeroline_bw -> Draw("same");
	
	c2 -> Update();
	
	
	theApp -> Run(kTRUE);
	
	
	//Saving part for the calibration functions
	string ans("");
	while(true){
		cout << "\nDo you want to save the calibration fit results? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans==string("y")  || ans==string("N") || ans==string("n") || ans==string("N")) break;
	}
	
	if(ans==string("y")  || ans==string("Y")){
		
		TFile *calibFile = new TFile((archivedir+calibset+string("/calibration.root")).c_str(),"UPDATE");
		
		calibFile -> WriteTObject(calib_fcn_fw,0,"overwrite");
		calibFile -> WriteTObject(calib_fcn_bw,0,"overwrite");
		calibFile -> WriteTObject(gr_cal_res_fw,"gr_cal_res_fw","overwrite");
		calibFile -> WriteTObject(gr_cal_res_bw,"gr_cal_res_bw","overwrite");
		
		calibFile -> Close();
	}
	
	
	
	//Resolution fitting
	
	c3 -> cd();
	gr_resol_ADC -> Draw("AP");
	
	TF1* resol_ADC_func = new TF1("resol_ADC_func","TMath::Sqrt([0]*x*x + [1]*x + [2])", (Double_t)gr_resol_ADC->GetXaxis()->GetFirst(), (Double_t)gr_resol_ADC->GetXaxis()->GetLast() );
	//TF1* resol_ADC_func = new TF1("resol_ADC_func", resolFunc, (Double_t)gr_resol_ADC->GetXaxis()->GetFirst(), (Double_t)gr_resol_ADC->GetXaxis()->GetLast(), 3);
	resol_ADC_func -> SetRange(0., 20000);
	resol_ADC_func -> SetParameters(1e-8,2e-3,5);
	resol_ADC_func -> SetLineColor(2);
	resol_ADC_func -> SetLineWidth(1);
	
	gr_resol_ADC -> Fit(resol_ADC_func);
	
	c3 -> Update();
	
	theApp -> Run(kTRUE);
	
	ans = "";
	while(true){
		cout << "\nDo you want to save the resolution fit results? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans==string("y")  || ans==string("N") || ans==string("n") || ans==string("N")) break;
	}
	
	if(ans==string("y")  || ans==string("Y")){
		TFile *resolFile = new TFile((archivedir+calibset+string("/resolution.root")).c_str(),"UPDATE");
		resolFile -> WriteTObject(resol_ADC_func, "resol_ADC_func", "Overwrite");
		resolFile -> WriteTObject(gr_resol_ADC, 0, "Overwrite");
		resolFile -> Close();
	}
	
	
	return 0;
}

void usage(){
	
	cout << '\nUsage:' << endl;
	cout << 'GatorFinalCalib <-d dataset>\n' << endl;
	cout << 'Note: a file named "fittedlines.root" (generated with "GatorBATCalib" prgram) should be present in the directory "' << archivedir << '/dataset".\n' << endl;
	
	return;
}
