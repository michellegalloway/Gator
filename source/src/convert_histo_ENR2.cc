#include <cstdlib>

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>

using namespace std;


TH1F* convert_histo_ENR2(TH1F* hADC, const char* calibdir){
	//This function gives the energy histogram with a constant binning. It is less accurate than convert_histo_ENR(...) routine but it is useful when I want build up the resolution function in energy domain.
	//cout <<"\nAcquisition time: " << time << " sec" << endl;
	
	TF1* calib_fcn_fw = 0;
	
	Int_t i=0; //General iterator index
	
	string calibfilename = calibdir;
	calibfilename += "calibration.root"; //Hardcoded but this function would be used for stuff referring to this.
	
	cout << "convert_histo_ENR2: Opening <" << calibfilename.c_str() << ">" << endl;
	TFile* f = new TFile(calibfilename.c_str(),"read");
	
	calib_fcn_fw = (TF1*)f->Get("calib_fcn_fw");
	
	
	if (calib_fcn_fw == 0){
		cout << "Calibration not still done. Run ./Th_calib before!!!" << endl;
		exit(1);
	}
	
	f -> Close();
	
	delete f;
	cout << "Calibration data loaded\n" << endl;
	
	//Using this constructor for the new histogram
	TH1F* hENR = new TH1F("hENR","",27500,0,2750); //each bin is 0.1 keV wide
	
	//Filling the new histogram
	Double_t enr_cent, counts;
	for (i=hADC->GetXaxis()->GetFirst(); i<=hADC->GetXaxis()->GetLast(); i++){
		
		enr_cent = calib_fcn_fw->Eval(hADC->GetBinCenter(i));
		counts = hADC->GetBinContent(i);
		hENR -> Fill(enr_cent,counts);
		//hENR -> SetBinError(hENR->FindBin(calib_fcn_fw->Eval(hADC->GetBinCenter(i))),tmp_rate_err);
		
	}
	
	return hENR; //The histogram which comes out must be rescaled for bin width -> do it outside this function
}
