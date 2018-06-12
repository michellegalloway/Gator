#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "TH1D.h"

#include "GatorStructs.h"


using namespace std;

///////////Definitions of the functions for the plot//////////////
Double_t peakFitFunc(Double_t* x,Double_t* par){
	
	Double_t E,P,T,A,sigma,beta,S,C;
	
	E = x[0];
	
	P = par[0];
	A = par[1];
	T = par[2];
	sigma = par[3];
	beta = par[4];
	S = par[5];
	C = par[6];
	
	return A*( (1-T)*TMath::Gaus(E,P,sigma) + T*TMath::Exp(beta*(E-P))*TMath::Erfc((E-(P-beta*sigma*sigma))/(sqrt(2)*sigma)) ) + S*TMath::Erfc((E-P)/(sqrt(2)*sigma)) + C;
	
	//return A*( TMath::Gaus(E,P,sigma) + R * TMath::Exp(beta*(E-P))*TMath::Erfc((E-(P-beta*sigma*sigma))/(sqrt(2)*sigma)) ) + C;
	
}
//////////////////////////////////////////////////////////////////////////////////


//----------Initialization functions for parameters-----------//

void costInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	
	Double_t left_mean = 0.;
	Double_t right_mean = 0.;
	
	for(int bin=firstbin; bin<firstbin+10; bin++){
		left_mean += (histo->GetBinContent(bin))/10.;
	}
	for(int bin=lastbin-9; bin<=lastbin; bin++){
		right_mean += (histo->GetBinContent(bin))/10.;
	}
	
	line.cost = TMath::Min(left_mean,right_mean);
	
	return;
}
//////////////////////////////////////////////////////////////////////////////////////

void stepInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	
	Double_t left_mean = 0.;
	Double_t right_mean = 0.;
	
	for(int bin=firstbin; bin<firstbin+10; bin++){
		left_mean += (histo->GetBinContent(bin))/10.;
	}
	for(int bin=lastbin-9; bin<=lastbin; bin++){
		right_mean += (histo->GetBinContent(bin))/10.;
	}
	
	if((left_mean - right_mean)<=0){
		line.step = right_mean;
	}else{
		line.step = (left_mean - right_mean);
	}
	
	return;
}
/////////////////////////////////////////////////////////////////////////////////////

void amplInit(TH1D* histo, CalibLine& line){
	//This should always be run after the function "costInit"
	double max =0.;
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	for(int bin=firstbin; bin<=lastbin; bin++){
		if(max < histo->GetBinContent(bin)) max = histo->GetBinContent(bin);
	}
	
	line.ampl = max-line.cost;
	//line.tail = 0.1*line.ampl;
	
	return;
}
//////////////////////////////////////////////////////////////////////////////////////

void sigmaInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	double max =0.;
	int maxbin;
	
	for(int bin=firstbin; bin<=lastbin; bin++){
		if(max < histo->GetBinContent(bin)){
			max = histo->GetBinContent(bin);
			maxbin = bin;
		}
	}
	
	int bin = maxbin;
	while(histo->GetBinContent(bin)>max/2.){
		bin--;
	}
	double lowval = histo->GetBinCenter(bin);
	
	bin = maxbin;
	while(histo->GetBinContent(bin)>max/2.){
		bin++;
	}
	double upval = histo->GetBinCenter(bin);
		
	line.sigma = (upval - lowval)/2; //This is not the sigma but a good starting point for the fit
	
	return;
}
/////////////////////////////////////////////////////////////////////////////////////


Double_t resolFunc (Double_t* x,Double_t* par){
	
	Double_t E = x[0];
	
	Double_t A = par[0];
	Double_t B = par[1];
	Double_t C = par[2];
	
	return TMath::Sqrt(A*E*E + B*E + C);
}