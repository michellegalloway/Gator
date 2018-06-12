#include <cstdlib>

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>

using namespace std;


TH1F* convert_histo_ENR(TH1F* hADC, const char* calibdir){
	//This function converts the histogram from ADC domain to KeV domain. But should be used only with samples acquisitions and backgrounds acquisitions, that refers to the current calibration, stored in the file calibration.root
	
	//cout <<"\nAcquisition time: " << time << " sec" << endl;
	
	string calibfilename = calibdir + string("calibration.root");
	
	Int_t i=0; //General iterator index
	
	//string calibfilename = calibdir;
	//calibfilename += "~/PhD/Gator/data_analysis/calibrations/calibration.root"; //Hardcoded but this function would be used for stuff referring to this.
	
	TFile* f = new TFile(calibfilename.c_str(),"read");
	
	TF1* calib_fcn_fw = (TF1*)f->Get("calib_fcn_fw");
	
	f -> Close();
	
	delete f;
	
	
	Double_t* xbinsENR = new Double_t[hADC->GetNbinsX()+1]; //Array of low energy bin edges (plus one more required by the contructor)
	
	for (i=0 ; i<hADC->GetNbinsX(); i++){
		xbinsENR[i] = calib_fcn_fw->Eval(hADC->GetXaxis()->GetBinLowEdge(i));
	}
	
	xbinsENR[hADC->GetNbinsX()] = calib_fcn_fw->Eval(hADC->GetXaxis()->GetBinUpEdge(i));// The upper edge of the last energy bin
	
	//Using this constructor for the new histogram
	TH1F* hENR = new TH1F("hENR","hENR",hADC->GetNbinsX(),xbinsENR);
	
	//Filling the new histogram
	for (i=hADC->GetXaxis()->GetFirst(); i<=hADC->GetXaxis()->GetLast(); i++){
				
		hENR -> Fill(calib_fcn_fw->Eval(hADC->GetBinCenter(i)),hADC->GetBinContent(i));
		//hENR -> SetBinError(hENR->FindBin(calib_fcn_fw->Eval(hADC->GetBinCenter(i))),tmp_rate_err);
		
	}
	
	delete xbinsENR;
	
	return hENR;
}
