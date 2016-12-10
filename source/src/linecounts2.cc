#include <TROOT.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLine.h>

#include <fstream>
#include <iostream>
#include <vector>


using namespace std;


void linecounts2(const TH1D* histoADC, //The units must be only in counts (not resdcaling)
				TCanvas* workcanv,
				const Double_t center,
				const Double_t sigma,
				const Double_t aqtime,
				Double_t& linecnts,
				Double_t& linecnts_err,
				Double_t& compcnts,
				Double_t& compcnts_err,
				Double_t& peakcnts,
				Double_t& peakcnts_err,
				Double_t& bgcnts, //This are "Line BG counts" (for peaked background)
				Double_t& bgcnts_err, //This are "Line BG counts" (for peaked background)
				Bool_t& peakedBG, //This is used outside if a peaked background is detected 
				Double_t& det_lim,
				const TH1F* bg_histoADC,//This is by default 0 ==> means that we are working with the BG only.
				Double_t bg_aqtime,//This is by default 0
				Double_t WidthScale){//This is by default 3.0
	
	//NOTE: here the bg_histoADC must enter in the routine with a properly rescaled for different time acquisitions.
	
	compcnts = 0;
	peakcnts = 0;
	bgcnts = 0;
	peakedBG = false;
	
	//The following six variables are for the 3 regions in ADC domain.
	Double_t line_low = center - WidthScale*sigma;
	Double_t line_up = center + WidthScale*sigma;
	Double_t line_intrv = line_up - line_low;
	Double_t compt_low = center - 2*WidthScale*sigma;
	Double_t compt_up = center + 2*WidthScale*sigma;
	Double_t compt_intrv = (line_low-compt_low) + (compt_up-line_up);
	
	
	//All the limits in general finish not on a bin border ==> 4 different bin fraction must be taken.
	//The fraction is relative to the lower section cut from the 4 lines
	Int_t bin1 = histoADC -> GetXaxis() -> FindBin(compt_low);
	Int_t bin2 = histoADC -> GetXaxis() -> FindBin(line_low);
	Int_t bin3 = histoADC -> GetXaxis() -> FindBin(line_up);
	Int_t bin4 = histoADC -> GetXaxis() -> FindBin(compt_up);
	
	
	Double_t frac1 = (compt_low-histoADC->GetXaxis()->GetBinLowEdge(bin1))/(histoADC->GetXaxis()->GetBinUpEdge(bin1)-histoADC->GetXaxis()->GetBinLowEdge(bin1));
	
	Double_t frac2 = (line_low-histoADC->GetXaxis()->GetBinLowEdge(bin2))/(histoADC->GetXaxis()->GetBinUpEdge(bin2)-histoADC->GetXaxis()->GetBinLowEdge(bin2));
	
	Double_t frac3 = (line_up-histoADC->GetXaxis()->GetBinLowEdge(bin3))/(histoADC->GetXaxis()->GetBinUpEdge(bin3)-histoADC->GetXaxis()->GetBinLowEdge(bin3));
	
	Double_t frac4 = (compt_up-histoADC->GetXaxis()->GetBinLowEdge(bin4))/(histoADC->GetXaxis()->GetBinUpEdge(bin4)-histoADC->GetXaxis()->GetBinLowEdge(bin4));
	
	
	compcnts += (1-frac1)*histoADC->GetBinContent(bin1);
	for (Int_t i=bin1+1; i<bin2; i++){
		compcnts += histoADC->GetBinContent(i);
	}
	compcnts += frac2*histoADC->GetBinContent(bin2);
	
	compcnts += (1-frac3)*histoADC->GetBinContent(bin3);
	for (Int_t i=bin3+1; i<bin4; i++){
		compcnts += histoADC->GetBinContent(i);
	}
	compcnts += frac4*histoADC->GetBinContent(bin4);
	
	peakcnts += (1-frac2)*histoADC->GetBinContent(bin2);
	for (Int_t i=bin2+1; i<bin3; i++){
		peakcnts += histoADC->GetBinContent(i);
	}
	peakcnts += frac3*histoADC->GetBinContent(bin3);
	
	//Up to here I have the gross peak counts in the variable "peakcnts" and the total Compton counts in the variable "peakcnts".
	
	
	//The following 3 variables will be used only if I am working with 2 histograms. Otherwise they remain unused.
	Double_t peakbg_detlim; //Threshold to see if a peak in the BG is detected
	Double_t bg_compcnts = 0; //Defined anyway but is not always used
	Double_t bg_compcnts_err; //Defined anyway but is not always used
	if(bg_histoADC){//For the BG counts reuse the code above and the same variables
		
		bin1 = bg_histoADC -> GetXaxis() -> FindBin(compt_low);
		bin2 = bg_histoADC -> GetXaxis() -> FindBin(line_low);
		bin3 = bg_histoADC -> GetXaxis() -> FindBin(line_up);
		bin4 = bg_histoADC -> GetXaxis() -> FindBin(compt_up);
		
		frac1 = (compt_low - bg_histoADC->GetXaxis()->GetBinLowEdge(bin1))/(bg_histoADC->GetXaxis()->GetBinUpEdge(bin1)-bg_histoADC->GetXaxis()->GetBinLowEdge(bin1));
		
		frac2 = (line_low - bg_histoADC->GetXaxis()->GetBinLowEdge(bin2)) / (bg_histoADC->GetXaxis()->GetBinUpEdge(bin2) - bg_histoADC->GetXaxis()->GetBinLowEdge(bin2));
	
		frac3 = (line_up - bg_histoADC->GetXaxis()->GetBinLowEdge(bin3)) / (bg_histoADC->GetXaxis()->GetBinUpEdge(bin3) - bg_histoADC->GetXaxis()->GetBinLowEdge(bin3));
		
		frac4 = (compt_up - bg_histoADC->GetXaxis()->GetBinLowEdge(bin4))/(bg_histoADC->GetXaxis()->GetBinUpEdge(bin4) - bg_histoADC->GetXaxis()->GetBinLowEdge(bin4));
	
		bgcnts += (1-frac2)*bg_histoADC->GetBinContent(bin2);
		for (Int_t i=bin2+1; i<bin3; i++){
			bgcnts += bg_histoADC->GetBinContent(i);
		}
		bgcnts += frac3*bg_histoADC->GetBinContent(bin3);
		//Up to here I have the gross BG peak counts in the variable bgcnts -> when transformed in net BG line counts the variable name will not be modified
		
		
		bg_compcnts += (1-frac1)*bg_histoADC->GetBinContent(bin1);
		for (Int_t i=bin1+1; i<bin2; i++){
			bg_compcnts += bg_histoADC->GetBinContent(i);
		}
		bg_compcnts += frac2*bg_histoADC->GetBinContent(bin2);
		
		bg_compcnts += (1-frac3)*bg_histoADC->GetBinContent(bin3);
		for (Int_t i=bin3+1; i<bin4; i++){
			bg_compcnts += bg_histoADC->GetBinContent(i);
		}
		bg_compcnts += frac4*bg_histoADC->GetBinContent(bin4);
		bg_compcnts_err = TMath::Sqrt(bg_compcnts);
		//Up to here I have the BG Compton counts in the variable "bg_compcnts"
		
		//Now for the BG net counts I subtract from "bgcnts" the estimated Compton and put the result in "bgcnts". But before I do the error estimation!!!
		bgcnts_err = TMath::Sqrt(bgcnts + bg_compcnts);
		bgcnts -= bg_compcnts; //From this point the variable "bgcnts" describes BG net counts
		
		peakbg_detlim = 2.86 + 4.87*TMath::Sqrt(1.36 + (line_intrv/compt_intrv)*bg_compcnts); //From the Gator paper;
		
		if(bgcnts>peakbg_detlim) peakedBG = true;
		
		bgcnts *= (aqtime/bg_aqtime); //Rescaled for sample acquisition time before using it
		bgcnts_err *= (aqtime/bg_aqtime); //Rescaled for sample acquisition time before using it
	}
	
	
	cout << "\nadc center = " << center << endl;
	cout << "adc sigma = " << sigma << endl;
	cout << "comptLow start = " << compt_low << endl;
	cout << "lineLow start = " << line_low << endl;
	cout << "lineUp stop = " << line_up << endl;
	cout << "comptUp stop = " << compt_up << endl;
	cout << "line_nbins = " << line_intrv << endl;
	cout << "compt_mbins = " << compt_intrv << endl;
	cout << "peakcnts = " << peakcnts << endl;
	cout << "peakcnts err = " << sqrt(peakcnts) << endl;
	cout << "raw comptcnts = " << compcnts << endl;
	cout << "raw comptcnts err = " << sqrt(compcnts) << endl;
	
	peakcnts_err = TMath::Sqrt(peakcnts);
	Double_t est_compcnts = (line_intrv/compt_intrv)*compcnts;
	Double_t est_compcnts_err = (line_intrv/compt_intrv)*sqrt(compcnts);
	
	cout << "est comptcnts = " << est_compcnts << endl;
	cout << "est comptcnts err = " << est_compcnts_err << endl;
	
	compcnts = est_compcnts;
	compcnts_err = est_compcnts_err;
	
	if (bg_histoADC==0){
		//We are working with a general BG acquired with the "signal" and it is used only to estimate a detection limit
		bgcnts = 0;
		bgcnts_err = 0;
	
		//In this situation the peak-counts are the BG itself
		det_lim = 2.86 + 4.87*TMath::Sqrt(1.36 + peakcnts); //From the Gator paper
	
		linecnts = peakcnts - est_compcnts;
		linecnts_err = TMath::Sqrt(peakcnts+(line_intrv/compt_intrv)*(line_intrv/compt_intrv)*compcnts); 
	
	
		cout << "linecnts = " << linecnts << endl;
		cout << "linecnts err = " << linecnts_err << endl;
	} else {
	//Here I have also the acquired "blanck" BG (properly rescaled)
		
		
		if(peakedBG){
			det_lim = 2.86 + 4.87*TMath::Sqrt(1.36 + (line_intrv/compt_intrv)*compcnts + bgcnts); //From the Gator paper
			linecnts = peakcnts - est_compcnts - bgcnts;
			linecnts_err = TMath::Sqrt(peakcnts + (line_intrv/compt_intrv)*(line_intrv/compt_intrv)*compcnts + bgcnts_err*bgcnts_err);
		}else{
			det_lim = 2.86 + 4.87*TMath::Sqrt(1.36 + (line_intrv/compt_intrv)*compcnts); //From the Gator paper
			linecnts = peakcnts - est_compcnts;
			linecnts_err = TMath::Sqrt(peakcnts + (line_intrv/compt_intrv)*(line_intrv/compt_intrv)*compcnts);
		}
		
		cout << "linecnts = " << linecnts << endl;
		cout << "linecnts err = " << linecnts_err << endl;
		
	}
	cout << "detect limit = " << det_lim << endl;
	
	if(bg_histoADC){
		//Output only if we deal with BG also!!!
		cout << "__Background section__" << endl;
		cout << "bg comptcnts (resc) = " << bg_compcnts*(aqtime/bg_aqtime) << endl;
		cout << "bg comptcnts err (resc) = " << bg_compcnts_err*(aqtime/bg_aqtime) << endl;
		cout << "bg net counts (resc) = " << bgcnts << endl;
		cout << "bg net counts err (resc) = " << bgcnts_err << endl;
		if(peakedBG) cout << "Peaked background!!!" << endl;
		
	}
	cout << "\n\n\n" << endl;
	
	
	TLine* tlinecenter = new TLine(center,workcanv->cd()->GetUymin(),center,workcanv->cd()->GetUymax());
	tlinecenter -> SetLineStyle(3);
	tlinecenter -> SetLineColor(kGreen+1);
	tlinecenter -> Draw("SAME");
	
	TLine* tline1 = new TLine(line_low,workcanv->cd()->GetUymin(),line_low,workcanv->cd()->GetUymax());
	tline1 -> SetLineStyle(3);
	tline1 -> SetLineColor(kRed);
	tline1 -> Draw("SAME");
	
	TLine* tline2 = new TLine(line_up,workcanv->cd()->GetUymin(),line_up,workcanv->cd()->GetUymax());
	tline2 -> SetLineStyle(3);
	tline2 -> SetLineColor(kRed);
	tline2 -> Draw("SAME");
	
	TLine* tline3 = new TLine(compt_low,workcanv->cd()->GetUymin(),compt_low,workcanv->cd()->GetUymax());
	tline3 -> SetLineStyle(3);
	tline3 -> SetLineColor(kBlue);
	tline3 -> Draw("SAME");
	
	TLine* tline4 = new TLine(compt_up,workcanv->cd()->GetUymin(),compt_up,workcanv->cd()->GetUymax());
	tline4 -> SetLineStyle(3);
	tline4 -> SetLineColor(kBlue);
	tline4 -> Draw("SAME");
	
	
	
	return;
}
