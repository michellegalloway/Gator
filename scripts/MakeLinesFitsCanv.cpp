#ifndef SCRIPT
#define SCRIPT

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <unistd.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
//#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
//#include "TCut.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
//#include "TIterator.h"
//#include "TList.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TApplication.h"
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>
//#include <TLatex.h>
//#include <TTimeStamp.h>
//#include <TLine.h>
#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
#include <TFitResult.h>

#include "../source/include/GatorStructs.h"
#include "../source/include/GatorCalibFitters.h"

using namespace std;


namespace Gator
{
	string CalibArchiveDir = "/Users/francesco/PhD/Gator/calibrations/archive/";
	
	vector<TCanvas*> canvases;
	vector<TH1D*> histograms;
	vector<TF1*> fitfunctions;
	vector<TGraphErrors*> residuals;
	
	TF1 *CalibFit=NULL, *ResolFit=NULL;
	
	void SpectralLinesFits(string date);
	
	void EnergyScaleFit(string date, bool resolfit=false);
	
}


void Gator::SpectralLinesFits(string date)
{
	
	TGaxis::SetMaxDigits(5);
	//gStyle->SetOptStat(0);
	gROOT->ForceStyle();
	
	chdir( (CalibArchiveDir+date).c_str() );
	
	int width=600, height=1000;
	double padsYborder=0.3, padLeftMargin=0.10, padRightMargin=0.05, labelsize1=0.03, titlesize1=0.04;
	
	string calibfilename = "fittedlines.root";
	
	TFile* calibfile = new TFile( calibfilename.c_str(), "read" );
	
	TTree* linestree = (TTree*)calibfile->Get("linestree");
	
	int nLines = linestree->GetEntries();
	
	
	TH1D* histo = NULL;
	//TF1* fit = NULL;
	string *mass= NULL, *element= NULL;
	double energy, mean, ampl, tail, sigma, beta, step, cost;
	
	
	linestree->SetBranchAddress("massNum", &mass);
	linestree->SetBranchAddress("element", &element);
	
	linestree->SetBranchAddress("litEn", &energy);
	
	//Parameters from fit results
	linestree->SetBranchAddress("mean",&mean);
	linestree->SetBranchAddress("ampl",&ampl);
	linestree->SetBranchAddress("tail",&tail);
	linestree->SetBranchAddress("sigma",&sigma);
	linestree->SetBranchAddress("beta",&beta);
	linestree->SetBranchAddress("step",&step);
	linestree->SetBranchAddress("cost",&cost);
	
	linestree->SetBranchAddress("histo", &histo);
	
	
	for(int iLine=0; iLine<nLines; iLine++){
		linestree->GetEntry(iLine);
		
		stringstream sstmp;
		
		sstmp.str(""); sstmp << (*element) << (*mass) << "_" << ( (int)(energy+0.5) ) << "keV_canv";
		TCanvas *canv = new TCanvas(sstmp.str().c_str(), sstmp.str().c_str(), width, height);
		canv->cd();
		
		
		sstmp.str(""); sstmp << (*element) << (*mass) << "_" << ( (int)(energy+0.5) ) << "keV_pad_" << 1;
		TPad *pad1 = new TPad(sstmp.str().c_str(), "", 0.0, padsYborder, 1.0, 1.0);
		pad1->SetLogy();
		pad1->SetTopMargin(0.01);
		pad1->SetBottomMargin(0);
		pad1->SetLeftMargin(padLeftMargin);
		pad1->SetRightMargin(padRightMargin);
		canv->cd(0);
		pad1->Draw();
		pad1->cd();
		
		sstmp << (*element) << (*mass) << "_" << ( (int)(energy+0.5) ) << "keV_histo";
		TH1D *h = (TH1D*)histo->Clone(sstmp.str().c_str());
		
		int nBins = h->GetNbinsX();
		double xMin = h->GetBinLowEdge(1);
		double xMax = h->GetBinLowEdge(nBins+1);
		TF1 *fit = new TF1("fit", peakFitFuncB, xMin, xMax, 7);
		fit->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
		fit->SetParameters(mean, ampl, tail, sigma, beta, step, cost);
		
		h->GetXaxis()->SetTickLength(0);
		h->GetXaxis()->SetTitleSize(titlesize1);
		h->GetXaxis()->SetTitleFont(132);
		h->GetXaxis()->SetLabelFont(132);
		h->GetXaxis()->SetLabelSize(labelsize1);
		
		h->GetYaxis()->SetTitleSize(titlesize1);
		h->GetYaxis()->SetTitleOffset(1.10);
		h->GetYaxis()->SetTitleFont(132);
		h->GetYaxis()->SetLabelFont(132);
		h->GetYaxis()->SetLabelSize(labelsize1);
		h->GetYaxis()->CenterTitle();
		
		h->Draw();
		fit->Draw("same");
		
		
		sstmp.str(""); sstmp << (*element) << (*mass) << "_" << ( (int)(energy+0.5) ) << "keV_pad_" << 2;
		TPad *pad2 = new TPad(sstmp.str().c_str(), "", 0.0, 0.0, 1.0, padsYborder);
		pad2->SetGridy();
		pad2->SetTopMargin(0);
		pad2->SetLeftMargin(padLeftMargin);
		pad2->SetBottomMargin(padLeftMargin*width/(height*padsYborder));
		pad2->SetRightMargin(padRightMargin);
		canv->cd(0);
		pad2->Draw();
		pad2->cd();
		
		TGraphErrors *gr = new TGraphErrors();
		gr->SetMarkerStyle(20);
		gr->SetMarkerSize(0.5);
		
		double yMax=0, yMin=0;
		for(int iBin=1; iBin<=nBins; iBin++){
			double x, ex, yex, y, ey;
			
			x = h->GetBinCenter(iBin);
			double xLow = h->GetBinLowEdge(iBin);
			double xUp = h->GetBinLowEdge(iBin+1);
			//ex = (xUp-xLow)/2;
			ex = 0;
			
			yex = (xUp-xLow)*(fit->Eval(xLow) + fit->Eval(xUp))/2;
			if(yex>0){
				y = (h->GetBinContent(iBin) - yex)/sqrt(yex);
				//ey = h->GetBinContent(iBin)/sqrt(yex);
				ey = 0;
				gr->SetPoint(iBin-1, x, y);
				gr->SetPointError(iBin-1, ex, ey);
				if(yMax<(y+ey)) yMax = (y+ey);
				if(yMin>(y-ey)) yMin = (y-ey);
			}
		}
		
		double yAbsMax = TMath::Max(yMax, TMath::Abs(yMin) );
		
		sstmp.str(""); sstmp << (*element) << (*mass) << "_" << ( (int)(energy+0.5) ) << "keV_frame";
		TH2D *frame = new TH2D(sstmp.str().c_str(), ";Channel;(y_{obs} - y_{exp}) / #sqrt{y_{exp}}", 100, xMin, xMax, 100, -1.1*yAbsMax, 1.1*yAbsMax);
		
		//This section is to make the same font heights in the two different pads
		double fs1 = h->GetXaxis()->GetTitleSize();
		double pw1 = pad1->XtoPixel(pad1->GetX2());
		double ph1 = pad1->YtoPixel(pad1->GetY1());
		double pw2 = pad2->XtoPixel(pad2->GetX2());
		double ph2 = pad2->YtoPixel(pad2->GetY1());
		double titlesize2 = titlesize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
		double labelsize2 = labelsize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
		
		frame->GetXaxis()->SetTickLength();
		frame->GetXaxis()->SetTitleSize(titlesize2);
		frame->GetXaxis()->SetTitleFont(132);
		frame->GetXaxis()->SetLabelFont(132);
		frame->GetXaxis()->SetLabelSize(labelsize2);
		frame->GetXaxis()->CenterTitle();
		
		frame->GetYaxis()->SetTickLength();
		frame->GetYaxis()->CenterTitle();
		frame->GetYaxis()->SetTitleOffset(0.4);
		frame->GetYaxis()->SetTitleSize(titlesize2);
		frame->GetYaxis()->SetTitleFont(132);
		frame->GetYaxis()->SetLabelFont(132);
		frame->GetYaxis()->SetLabelSize(labelsize2);
		frame->GetYaxis()->CenterTitle();
		
		frame->Draw();
		
		
		gr->SetMarkerColor(kBlack);
		gr->SetLineColor(kBlack);
		gr->SetLineWidth(1);
		gr->Draw("P");
		
		TF1 *zeroline = new TF1("zeroline", "0*x", xMin, xMax);
		zeroline->SetLineColor(kRed);
		zeroline->Draw("same");
		
	}
}


void Gator::EnergyScaleFit(string date, bool resolfit)
{
	TGaxis::SetMaxDigits(5);
	//gStyle->SetOptStat(0);
	gROOT->ForceStyle();
	
	TVirtualFitter::SetDefaultFitter("Minuit2");
	
	chdir( (CalibArchiveDir+date).c_str() );
	
	int width=800, height=800;
	double padsYborder=0.33, padLeftMargin=0.10, padRightMargin=0.05, labelsize1=0.03, titlesize1=0.04;
	
	string calibfilename = "fittedlines.root";
	
	TFile* calibfile = new TFile(calibfilename.c_str(), "read");
	
	TTree* linestree = (TTree*)calibfile->Get("linestree");
	
	int nLines = linestree->GetEntries();
	
	
	double energy, energy_err, mean, mean_err, sigma, sigma_err;
	
	
	linestree->SetBranchAddress("litEn", &energy);
	linestree->SetBranchAddress("litEn_err", &energy_err);
	linestree->SetBranchAddress("mean", &mean);
	linestree->SetBranchAddress("mean_err", &mean_err);
	linestree->SetBranchAddress("sigma", &sigma);
	linestree->SetBranchAddress("sigma_err", &sigma_err);
	
	
	TGraphErrors* grCal = new TGraphErrors();
	grCal->SetMarkerStyle(7);
	grCal->SetMarkerSize(2);
	
	double xMin, xMax, yMin, yMax;
	
	for(int iLine=0; iLine<nLines; iLine++){
		
		linestree->GetEntry(iLine);
		
		grCal->SetPoint(iLine, mean, energy);
		grCal->SetPointError(iLine, mean_err, energy_err);
		
		if(iLine==0){
			xMin = mean-mean_err;
			xMax = mean+mean_err;
			yMin = energy-energy_err;
			yMax = energy+energy_err;
		}else{
			if(xMin>(mean-mean_err)) xMin = mean-mean_err;
			if(xMax<(mean+mean_err)) xMax = mean+mean_err;
			if(yMin>(energy-energy_err)) yMin = energy-energy_err;
			if(yMax<(energy+energy_err)) yMax = energy+energy_err;
		}
	}
	
	double Dx = xMax-xMin;
	double Dy = yMax-yMin;
	
	TCanvas *canv1 = new TCanvas("canv1","Energy calibration", width, height);
	
	canv1->cd(0);
	TPad *pad1_1 = new TPad("pad1_1", "", 0.0, padsYborder, 1.0, 1.0);
	pad1_1->SetTopMargin(0.01);
	pad1_1->SetBottomMargin(0);
	pad1_1->SetLeftMargin(padLeftMargin);
	pad1_1->SetRightMargin(padRightMargin);
	pad1_1->Draw();
	pad1_1->cd();
	
	TH1F *frame1_1 = pad1_1->DrawFrame(xMin-(Dx/20), yMin-(Dy/20), xMax+(Dx/20), yMax+(Dy/20), ";Channel;Energy [keV]" );
	
	frame1_1->GetXaxis()->CenterTitle();
	frame1_1->GetXaxis()->SetTickLength(0);
	frame1_1->GetXaxis()->SetTitleSize(titlesize1);
	frame1_1->GetXaxis()->SetTitleFont(132);
	frame1_1->GetXaxis()->SetLabelFont(132);
	frame1_1->GetXaxis()->SetLabelSize(labelsize1);
	
	frame1_1->GetYaxis()->CenterTitle();
	frame1_1->GetYaxis()->SetTitleSize(titlesize1);
	frame1_1->GetYaxis()->SetTitleOffset(1.05);
	frame1_1->GetYaxis()->SetTitleFont(132);
	frame1_1->GetYaxis()->SetLabelFont(132);
	frame1_1->GetYaxis()->SetLabelSize(labelsize1);
	frame1_1->GetYaxis()->SetDecimals();
	
	grCal->Draw("P");
	
	if(CalibFit){
		delete CalibFit;
		CalibFit = NULL;
	}
	CalibFit = new TF1("CalibFit", "pol2", xMin, xMax);
	TFitResultPtr result = grCal->Fit(CalibFit,"MVS");
	pad1_1->Update();
	
	
	canv1->cd(0);
	TPad *pad1_2 = new TPad("pad1_2", "", 0.0, 0.0, 1.0, padsYborder);
	pad1_2->SetGridy();
	pad1_2->SetTopMargin(0);
	pad1_2->SetLeftMargin(padLeftMargin);
	pad1_2->SetBottomMargin(padLeftMargin*width/(height*padsYborder));
	pad1_2->SetRightMargin(padRightMargin);
	pad1_2->Draw();
	pad1_2->cd();
	
	TGraphErrors* grCalResid = new TGraphErrors();
	grCalResid->GetXaxis()->SetTitle("Channel"); grCalResid->GetXaxis()->CenterTitle();
	grCalResid->GetYaxis()->SetTitle("Residuals [KeV]"); grCalResid->GetYaxis()->CenterTitle();
	grCalResid->SetMarkerStyle(7);
	grCalResid->SetMarkerSize(2);
	
	
	for(int iLine=0; iLine<nLines; iLine++){
		linestree->GetEntry(iLine);
		
		double fitVal = CalibFit->Eval(mean);
		grCalResid->SetPoint(iLine, mean, energy-fitVal);
		//grCalRes->SetPointError(iLine, mean_err, energy_err);
		
		if(iLine==0){
			yMin = energy-fitVal;
			yMax = energy-fitVal;
		}else{
			if(yMin>(energy-fitVal)) yMin = energy-fitVal;
			if(yMax<(energy-fitVal)) yMax = energy-fitVal;
		}
		
	}
	
	Dy = yMax-yMin;
	
	
	TH1F *frame1_2 = pad1_2->DrawFrame(xMin-(Dx/20), yMin-(Dy/10), xMax+(Dx/20), yMax+(Dy/10), ";Channel;Residuals [keV]" );
	
	{
	//This section is to make the same font heights in the two different pads
	//double fs1 = frame1_2->GetXaxis()->GetTitleSize();
	double pw1 = pad1_1->XtoPixel(pad1_1->GetX2());
	double ph1 = pad1_1->YtoPixel(pad1_1->GetY1());
	double pw2 = pad1_2->XtoPixel(pad1_2->GetX2());
	double ph2 = pad1_2->YtoPixel(pad1_2->GetY1());
	double titlesize2 = titlesize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
	double labelsize2 = labelsize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
	
	frame1_2->GetXaxis()->CenterTitle();
	frame1_2->GetXaxis()->SetTickLength(0);
	frame1_2->GetXaxis()->SetTitleSize(titlesize2);
	frame1_2->GetXaxis()->SetTitleFont(132);
	frame1_2->GetXaxis()->SetLabelFont(132);
	frame1_2->GetXaxis()->SetLabelSize(labelsize2);
	
	frame1_2->GetYaxis()->CenterTitle();
	frame1_2->GetYaxis()->SetTitleSize(titlesize2);
	frame1_2->GetYaxis()->SetTitleOffset(0.5);
	frame1_2->GetYaxis()->SetTitleFont(132);
	frame1_2->GetYaxis()->SetLabelFont(132);
	frame1_2->GetYaxis()->SetLabelSize(labelsize2);
	frame1_2->GetYaxis()->SetDecimals();
	}
	grCalResid->Draw("P");
	
	TF1 *zeroline = new TF1("zeroline", "0*x", 0, 1e6);
	zeroline->SetLineColor(kRed);
	zeroline->Draw("same");
	
	
	
	if(!resolfit) return;
	
	cout << "\n\n\nENERGY RESOLUTION\n\n\n" << endl;
	
	TCanvas *canv2 = new TCanvas("canv2","Energy resolution", width, height);
	canv2->cd(0);
	TPad *pad2_1 = new TPad("pad2_1", "", 0.0, padsYborder, 1.0, 1.0);
	pad2_1->SetTopMargin(0.01);
	pad2_1->SetBottomMargin(0);
	pad2_1->SetLeftMargin(padLeftMargin);
	pad2_1->SetRightMargin(padRightMargin);
	pad2_1->Draw();
	pad2_1->cd();
	
	
	//Fit parameters
	double a = CalibFit->GetParameter(2);
	double b = CalibFit->GetParameter(1);
	double c = CalibFit->GetParameter(0);
	
	//Parameter errors
	double a_err = CalibFit->GetParError(2);
	double b_err = CalibFit->GetParError(1);
	double c_err = CalibFit->GetParError(0);
	
	//Covariance elements matrix 
	double Cov_ab = result->CovMatrix(1,2);
	double Cov_bc = result->CovMatrix(1,0);
	double Cov_ac = result->CovMatrix(0,2);
	
	TGraphErrors *grResol = new TGraphErrors();
	
	for(int iLine=0; iLine<nLines; iLine++){
		linestree->GetEntry(iLine);
		
		//Clculate the sigma in energy units
		double sE = sigma*sqrt( 2*pow(a*sigma,2) + pow(2*a*a*mean,2) + 4*a*b*mean + pow(b,2) );
		
		//Partial derivatives for error propagation
		double dsE_dmean = (pow(sigma,2)/sE)*( 4*pow(a,2)*mean + 2*a*b );
		double dsE_dsigma = (pow(sigma,2)/sE)*( 4*pow(a,2)*sigma );
		double dsE_da = (pow(sigma,2)/sE)*( 2*a*pow(sigma,2) + 4*a*pow(mean,2) + 2*b*mean );
		double dsE_db = (pow(sigma,2)/sE)*( 2*a*mean + b );
		
		//Calculate the error on the sigma in energy units
		//double sE_err = sqrt( pow(dsE_dmean*mean_err,2) + pow(dsE_dsigma*sigma_err,2) + pow(dsE_da*a_err,2) + pow(dsE_db*b_err,2) + 2*dsE_da*dsE_db*Cov_ab ); //This formula underestimates the errors
		
		double sE_err = (2*a*mean + b)*sigma_err;//Based on the propagation of an interval with a relatively linaer function
		
		grResol->SetPoint(iLine, energy, sE);
		grResol->SetPointError(iLine, energy_err, sE_err);
		
		if(iLine==0){
			xMin = energy-energy_err;
			xMax = energy+energy_err;
			yMin = sE-sE_err;
			yMax = sE+sE_err;
		}else{
			if(xMin>(energy-energy_err)) xMin = energy-energy_err;
			if(xMax<(energy+energy_err)) xMax = energy+energy_err;
			if(yMin>(sE-sE_err)) yMin = sE-sE_err;
			if(yMax<(sE+sE_err)) yMax = sE+sE_err;
		}
		
		//cout << "\nEnrgy = " << energy << " keV" << endl;
		//cout << "Sigma(E) = " << sE << " keV" << endl;
		//cout << "D[Sigma(E)] " << sE_err << " keV" << endl;
	}
	
	Dx = xMax-xMin;
	Dy = yMax-yMin;
	
	TH1F *frame2_1 = pad2_1->DrawFrame(xMin-(Dx/20), yMin-(Dy/20), xMax+(Dx/20), yMax+(Dy/20), ";Energy [keV];Energy resolution [keV]" );
	
	frame2_1->GetXaxis()->CenterTitle();
	frame2_1->GetXaxis()->SetTickLength(0);
	frame2_1->GetXaxis()->SetTitleSize(titlesize1);
	frame2_1->GetXaxis()->SetTitleFont(132);
	frame2_1->GetXaxis()->SetLabelFont(132);
	frame2_1->GetXaxis()->SetLabelSize(labelsize1);
	
	frame2_1->GetYaxis()->CenterTitle();
	frame2_1->GetYaxis()->SetTitleSize(titlesize1);
	frame2_1->GetYaxis()->SetTitleOffset(1.05);
	frame2_1->GetYaxis()->SetTitleFont(132);
	frame2_1->GetYaxis()->SetLabelFont(132);
	frame2_1->GetYaxis()->SetLabelSize(labelsize1);
	frame2_1->GetYaxis()->SetDecimals();
	
	grResol->Draw("P");
	
	
	if(ResolFit){
		delete ResolFit;
		ResolFit = NULL;
	}
	ResolFit = new TF1("ResolFit", "sqrt([0] + [1]*x + [2]*pow(x,2))+[3]", xMin, xMax);
	ResolFit->SetParameters(1e-1, 3e-4, 1e-8, 0);
	ResolFit->FixParameter(0, 0.);
	//ResolFit->FixParameter(2, 0.);
	//ResolFit->FixParameter(3, 0.);
	TFitResultPtr result2 = grResol->Fit(ResolFit,"MVS");
	pad2_1->Update();
	
	
	
	
	canv2->cd(0);
	TPad *pad2_2 = new TPad("pad2_2", "", 0.0, 0.0, 1.0, padsYborder);
	pad2_2->SetGridy();
	pad2_2->SetTopMargin(0);
	pad2_2->SetLeftMargin(padLeftMargin);
	pad2_2->SetBottomMargin(padLeftMargin*width/(height*padsYborder));
	pad2_2->SetRightMargin(padRightMargin);
	pad2_2->Draw();
	pad2_2->cd();
	
	TGraphErrors* grResolResid = new TGraphErrors();
	grResolResid->GetXaxis()->SetTitle("Energy [keV]"); grResolResid->GetXaxis()->CenterTitle();
	grResolResid->GetYaxis()->SetTitle("Residuals [KeV]"); grResolResid->GetYaxis()->CenterTitle();
	grResolResid->SetMarkerStyle(7);
	grResolResid->SetMarkerSize(2);
	
	
	for(int iLine=0; iLine<nLines; iLine++){
		linestree->GetEntry(iLine);
		
		//Calculate the sigma in energy units
		double sE = sigma*sqrt( 2*pow(a*sigma,2) + pow(2*a*a*mean,2) + 4*a*b*mean + pow(b,2) );
		
		//Partial derivatives for error propagation
		double dsE_dmean = (pow(sigma,2)/sE)*( 4*pow(a,2)*mean + 2*a*b );
		double dsE_dsigma = (pow(sigma,2)/sE)*( 4*pow(a,2)*sigma );
		double dsE_da = (pow(sigma,2)/sE)*( 2*a*pow(sigma,2) + 4*a*pow(mean,2) + 2*b*mean );
		double dsE_db = (pow(sigma,2)/sE)*( 2*a*mean + b );
		
		//Calculate the error on the sigma in energy units
		//double sE_err = sqrt( pow(dsE_dmean*mean_err,2) + pow(dsE_dsigma*sigma_err,2) + pow(dsE_da*a_err,2) + pow(dsE_db*b_err,2) + 2*dsE_da*dsE_db*Cov_ab ); //This formula underestimates the errors
		
		double sE_err = (2*a*mean + b)*sigma_err;//Based on the propagation of an interval with a relatively linaer function
		
		
		double fitVal = ResolFit->Eval(energy);
		grResolResid->SetPoint(iLine, energy, sE-fitVal);
		//grCalRes->SetPointError(iLine, mean_err, energy_err);
		
		if(iLine==0){
			yMin = sE-fitVal;
			yMax = sE-fitVal;
		}else{
			if(yMin>(sE-fitVal)) yMin = sE-fitVal;
			if(yMax<(sE-fitVal)) yMax = sE-fitVal;
		}
		
	}
	
	Dy = yMax-yMin;
	
	
	TH1F *frame2_2 = pad2_2->DrawFrame(xMin-(Dx/20), yMin-(Dy/10), xMax+(Dx/20), yMax+(Dy/10), ";Energy [keV];Residuals [keV]" );
	
	{
		//This section is to make the same font heights in the two different pads
		double pw1 = pad2_1->XtoPixel(pad2_1->GetX2());
		double ph1 = pad2_1->YtoPixel(pad2_1->GetY1());
		double pw2 = pad2_2->XtoPixel(pad2_2->GetX2());
		double ph2 = pad2_2->YtoPixel(pad2_2->GetY1());
		double titlesize2 = titlesize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
		double labelsize2 = labelsize1*TMath::Min(pw1,ph1)/TMath::Min(pw2,ph2);
	
		frame2_2->GetXaxis()->CenterTitle();
		frame2_2->GetXaxis()->SetTickLength(0);
		frame2_2->GetXaxis()->SetTitleSize(titlesize2);
		frame2_2->GetXaxis()->SetTitleFont(132);
		frame2_2->GetXaxis()->SetLabelFont(132);
		frame2_2->GetXaxis()->SetLabelSize(labelsize2);
	
		frame2_2->GetYaxis()->CenterTitle();
		frame2_2->GetYaxis()->SetTitleSize(titlesize2);
		frame2_2->GetYaxis()->SetTitleOffset(0.5);
		frame2_2->GetYaxis()->SetTitleFont(132);
		frame2_2->GetYaxis()->SetLabelFont(132);
		frame2_2->GetYaxis()->SetLabelSize(labelsize2);
		frame2_2->GetYaxis()->SetDecimals();
	}
	
	grResolResid->Draw("P");
	zeroline->DrawCopy("same");
}

#endif


#include "../source/src/GatorCalibFitters.cc"

