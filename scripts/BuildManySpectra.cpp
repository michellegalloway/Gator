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
//#include <TH1D.h>
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
TH1D* loadSpe(const char* dir, Double_t& aqtime);
TH1D* convert_histo_ENR(TH1D* hADC, TF1* calibfunc);



namespace GatorSpectra
{
	class Spectrum{
	public:
		Spectrum()
		{
			name = "";
			RawSpec = NULL;
			EnSpect = NULL;
			EnSpecScaled = NULL;
			Calib = NULL;
			RawResol = NULL;
			EnResol = NULL;
			AqTime = 0;
			color = kBlack;
		};
		
		Spectrum(string _name, TH1D* spec, double _AqTime, Color_t _color=kBlack)
		{
			name = _name;
			RawSpec = spec;
			
			EnSpect = NULL;
			EnSpecScaled = NULL;
			Calib = NULL;
			RawResol = NULL;
			EnResol = NULL;
			
			AqTime = _AqTime;
			color = _color;
		};
		
		void SetRawSpectrum(string _name, TH1D* spec, double _AqTime, Color_t _color=kBlack);
		void LoadCalib(string calibdir);
		
		string GetName(){return name;};
		
		TH1D* GetHisto(){return EnSpecScaled;};
		TH1D* GetHistoEnNotScaled(){return EnSpect;};
		Color_t GetColor(){return color;};
		
		void SetCalib(TF1 *calib);
		void SetRawResol(TF1 *resol){;};
		void SetEnResol(TF1 *resol){;};
		
		//void SetRebin(int rebin){rebinning=rebin;};
		
		void MakeCalibSpectrum();
		
		void Draw(string opt="");
		
	private:
		string name;
		TH1D *RawSpec, *EnSpecScaled, *EnSpect;
		TF1 *Calib, *RawResol, *EnResol;
		double AqTime;
		Color_t color;
	};
}


namespace GatorSpectra
{
	
	string SamplesArchive = "/Users/francesco/PhD/Gator/analysis/samples/";
	string BkgdArchiveDir = "/Users/francesco/PhD/Gator/background/archive/";
	
	string CalibArchiveDir = "/Users/francesco/PhD/Gator/calibrations/archive/";
	
	
	
	
	
	void AddSpectrum(string name, string dir, string calib="", Color_t _color=kBlack);
	void AddCalibSourceSpectrum(string name, string dir, string calib="", Color_t _color=kBlack);
	void SetBkgdSpectrum(string dir, string calib="");
	
	void SetUpDrawingStyle();
	
	void MakeCalibSpectra();
	void MakeSteelSpectra();
	void MakeTitaniumSpectra();
	void MakePTFESpectra();
	void MakeCopperSpectra();
	void MakeCeramicStemSpectra();
	void MakePMTholderSpectra();
	void MakePMTsBatchSpectra();
	void MakePMTsFlangesSpectra();
	
	vector<Spectrum*> SpectraVec;
	Spectrum* BkgnSpect=NULL;
	
	
	TCanvas* c1=NULL;
	
	int rebinning = 16;
	
	double minEn = 50;
	double maxEn = 2700;
	
	bool makeupsettings=false;
}



void GatorSpectra::Spectrum::SetRawSpectrum(string _name, TH1D* spec, double _AqTime, Color_t _color)
{
	if(RawSpec) delete RawSpec;
	if(EnSpect) delete EnSpect;
	if(EnSpecScaled) delete EnSpecScaled;
	if(Calib) delete Calib;
	if(RawResol) delete RawResol;
	if(EnResol) delete EnResol;
	
	name = _name;
	RawSpec = spec;
	RawSpec->SetName( (name+string("Raw")).c_str() );
	AqTime = _AqTime;
}


void GatorSpectra::Spectrum::SetCalib(TF1 *calib)
{
	if(Calib) delete Calib;
	Calib = calib;
}


void GatorSpectra::Spectrum::MakeCalibSpectrum()
{
	if(!Calib) return;
	if(!RawSpec) return;
	
	if(EnSpecScaled) delete EnSpecScaled;
	if(EnSpect) delete EnSpect;
	
	
	EnSpect = convert_histo_ENR(RawSpec, Calib);
	EnSpect->SetName( (name+string("Energy")).c_str() );
	
	EnSpecScaled = (TH1D*)EnSpect->Clone( (name+string("Scaled")).c_str() );
	EnSpecScaled->Rebin( GatorSpectra::rebinning );
	EnSpecScaled->Scale( (24*3600.)/AqTime, "width" );
	EnSpecScaled->SetTitle(";Energy [keV];Differential rate [counts keV^{-1} day^{-1}]");
	EnSpecScaled->SetLineColor(color);
	
	EnSpecScaled->GetXaxis()->SetLabelFont(132); //Times New Roman
	EnSpecScaled->GetXaxis()->SetLabelSize(0.045);
	EnSpecScaled->GetXaxis()->SetTitleFont(132); //Times New Roman
	EnSpecScaled->GetXaxis()->SetTitleSize(0.045);
	EnSpecScaled->GetXaxis()->SetTitleOffset(1.1);
	
	EnSpecScaled->GetYaxis()->SetLabelFont(132); //Times New Roman
	EnSpecScaled->GetYaxis()->SetLabelSize(0.045);
	EnSpecScaled->GetYaxis()->SetTitleFont(132); //Times New Roman
	EnSpecScaled->GetYaxis()->SetTitleSize(0.045);
	EnSpecScaled->GetYaxis()->SetTitleOffset(0.55);
	
}


void GatorSpectra::Spectrum::LoadCalib(string calibdir)
{
	string filename = GatorSpectra::CalibArchiveDir + calibdir + string("/calibration.root");
	
	TFile *f = TFile::Open( filename.c_str(),"read");
	
	if(!f) return;
	
	TF1* calib = (TF1*)f->Get("calib_fcn_fw");
	if(!calib){
		cerr << "\nERROR --> GatorSpectra::Spectrum::LoadCalib(...): Cannot find the calibration function in file <" << f->GetName() << ">" << endl;
		delete f;
		return;
	}
	
	SetCalib( (TF1*)calib->Clone("Calib") );
	
	delete f;
	
}


void GatorSpectra::Spectrum::Draw(string opt)
{
	GatorSpectra::SetUpDrawingStyle();
	
	if(!EnSpecScaled) return;
	
	//EnSpecScaled->SetLineWidth(2);
	EnSpecScaled->GetXaxis() -> SetRangeUser( minEn, maxEn );
	EnSpecScaled->Draw( opt.c_str() );
	
}




void GatorSpectra::AddCalibSourceSpectrum(string name, string dir, string calib, Color_t _color)
{

	double aqtime;
	string spedir = CalibArchiveDir+dir;
	TH1D* hRaw = loadSpe(spedir.c_str(), aqtime);
	GatorSpectra::Spectrum *spec = new GatorSpectra::Spectrum(name, hRaw, aqtime, _color);
	spec->LoadCalib(calib);
	spec->MakeCalibSpectrum();
	
	
	unsigned nDs = SpectraVec.size();
	
	bool present = false;
	for(int iDs=0; iDs<nDs; iDs++){
		if(name == SpectraVec.at(iDs)->GetName()){
			present=true;
			delete SpectraVec.at(iDs);
			SpectraVec.at(iDs) = spec;
		}
	}
	
	if(!present) SpectraVec.push_back(spec);
	
	return;
}


void GatorSpectra::AddSpectrum(string name, string dir, string calib, Color_t _color)
{

	double aqtime;
	string spedir = SamplesArchive+dir;
	TH1D* hRaw = loadSPE(spedir.c_str(), aqtime);
	GatorSpectra::Spectrum *spec = new GatorSpectra::Spectrum(name, hRaw, aqtime, _color);
	spec->LoadCalib(calib);
	spec->MakeCalibSpectrum();
	
	
	unsigned nDs = SpectraVec.size();
	
	bool present = false;
	for(int iDs=0; iDs<nDs; iDs++){
		if(name == SpectraVec.at(iDs)->GetName()){
			present=true;
			delete SpectraVec.at(iDs);
			SpectraVec.at(iDs) = spec;
		}
	}
	
	if(!present) SpectraVec.push_back(spec);
	
	return;
}


void GatorSpectra::SetBkgdSpectrum(string dir, string calibdir)
{
	if(BkgnSpect) delete BkgnSpect;
	
	double aqtime;
	string spedir = BkgdArchiveDir+dir;
	TH1D* hRaw = loadSPE(spedir.c_str(), aqtime);
	BkgnSpect = new GatorSpectra::Spectrum("Background", hRaw, aqtime, kGray+2);
	BkgnSpect->LoadCalib(calibdir);
	BkgnSpect->MakeCalibSpectrum();
	BkgnSpect->GetHisto()->SetLineWidth(2);
	
}


void GatorSpectra::SetUpDrawingStyle()
{
	if(makeupsettings) return;
	
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
	
	makeupsettings=true;
}


//========================================//

void GatorSpectra::MakeCalibSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.62);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddCalibSourceSpectrum("^{228}Th source", "2013.07.01/228Th/", "2013.07.01", kBlue-2);
	GatorSpectra::AddCalibSourceSpectrum("^{60}Co source ", "2015.08.07/60Co/", "2015.08.07", kRed-2);
	
	unsigned nDs=SpectraVec.size();
	
	for(unsigned iDs=0; iDs<nDs; iDs++){
		
		//Scale by area the spectra
		TH1D *hist = SpectraVec.at(iDs)->GetHistoEnNotScaled();
		hist->SetLineColor(SpectraVec.at(iDs)->GetColor());
		hist->Rebin(rebinning);
		double integral = hist->Integral(minEn,maxEn);
		hist->SetTitle(";Energy [keV];");
		hist->Scale(1./integral);
		
		if(iDs==0){
			hist->Draw();
			leg->AddEntry(hist, (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			hist->Draw("same");
			leg->AddEntry(hist, (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	
	leg->Draw();
}


void GatorSpectra::MakeSteelSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.62);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("Sample 1", "NironitSS_2012/SPE/", "2012.10.04", kBlue-3);
	GatorSpectra::AddSpectrum("Sample 4", "25cmSteelTube/SPE/", "2013.12.23", kRed-3);
	GatorSpectra::AddSpectrum("Sample 5", "ssteel_pipes/SPE/", "2012.10.04", kGreen-2);
	GatorSpectra::AddSpectrum("Sample 7", "NironitSS_2013/SPE/", "2013.12.23", kOrange+2);
	GatorSpectra::AddSpectrum("Sample 17", "WTostoSS/SPE/", "2012.10.04", kCyan-3);
	GatorSpectra::AddSpectrum("Sample 20", "SS_Alca_09072013/SPE/", "2013.07.01", kViolet+1);
	
	GatorSpectra::SetBkgdSpectrum("2012/", "2012.06.19");
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			SpectraVec.at(iDs)->Draw();
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	
	leg->Draw();
	
}


void GatorSpectra::MakeTitaniumSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("Sample 1", "TiGr1_Uniti/SPE/", "2012.06.19", kBlue-3);
	GatorSpectra::AddSpectrum("Sample 2", "TiGr1_Thyssen/SPE/", "2012.06.19", kRed-3);
	
	GatorSpectra::SetBkgdSpectrum("2012/", "2012.06.19");
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			SpectraVec.at(iDs)->Draw();
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	
	leg->Draw();
	
}


void GatorSpectra::MakePTFESpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("PTFE complete run", "PTFE_40kg/SPE/", "2013.12.23", kRed-3);
	GatorSpectra::AddSpectrum("PTFE reduced run", "PTFE_40kg/SPEreduced/", "2013.12.23", kBlue-3);
	
	GatorSpectra::SetBkgdSpectrum("2012/", "2012.06.19");
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	/*
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	*/
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
}


void GatorSpectra::MakeCopperSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("Gold plated OFHC Cu", "AuPlatedCuRods/SPE/", "2013.12.23", kRed-3);
	
	GatorSpectra::SetBkgdSpectrum("2014/", "2013.12.23");
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	/*
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	*/
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
}


void GatorSpectra::MakeCeramicStemSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("Ceramic (Al_{2}O_{3})", "PMTs/R11410-21/BulkMaterials/ceramic/SPE/", "2012.10.04", kRed-3);
	
	GatorSpectra::SetBkgdSpectrum("2014/", "2013.12.23");
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	/*
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	*/
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
}


void GatorSpectra::MakePMTholderSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("PTFE holders 2013", "PTFE_holders/SPE/", "2013.07.01", kRed-3);
	GatorSpectra::AddSpectrum("PTFE holders 2014", "PTFE_holders/SPE2/", "2013.12.23", kBlue+2);
	
	GatorSpectra::SetBkgdSpectrum("2014/", "2013.12.23");
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	/*
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	*/
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
}


void GatorSpectra::MakePMTsBatchSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	
	
	GatorSpectra::AddSpectrum("PMTs batch 2", "PMTs/R11410-21/Batch2/SPE/", "2013.07.01", kRed-3);
	GatorSpectra::AddSpectrum("PMTs batch 8", "PMTs/R11410-21/Batch8/SPE/", "2013.12.23", kBlue+2);
	GatorSpectra::AddSpectrum("PMTs batch 20", "PMTs/R11410-21/Batch20/SPE/", "2014.10.10", kGreen-2);
	
	GatorSpectra::AddSpectrum("PTFE holders 2014", "PTFE_holders/SPE2/", "2013.12.23", kCyan-6);
	
	GatorSpectra::SetBkgdSpectrum("2014/", "2013.12.23");
	
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
	
	return;
}


void GatorSpectra::MakePMTsFlangesSpectra()
{
	if(c1) delete c1;
	
	c1 = new TCanvas("c1","",1200,500);
	c1->SetLogy();

	
	TLegend* leg = new TLegend(0.65, 0.85, 0.8, 0.60);
	leg->SetBorderSize(0);
	
	GatorSpectra::AddSpectrum("R11410-10 flanges", "PMTs/CeramicStems_R11410-10/SPE/", "2013.07.01", kRed-3);
	GatorSpectra::AddSpectrum("R11410-21 flanges", "PMTs/CeramicStems_R11410-21/SPE/", "2013.07.01", kBlue+2);
	
	GatorSpectra::SetBkgdSpectrum("2014/", "2013.12.23");
	
	if(BkgnSpect){
		BkgnSpect->Draw();
	}
	
	
	unsigned nDs=SpectraVec.size();
	for(unsigned iDs=0; iDs<nDs; iDs++){
		if(iDs==0){
			if(BkgnSpect){
				SpectraVec.at(iDs)->Draw("same");
			}else{
				SpectraVec.at(iDs)->Draw();
			}
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}else{
			SpectraVec.at(iDs)->Draw("same");
			leg->AddEntry(SpectraVec.at(iDs)->GetHisto(), (SpectraVec.at(iDs)->GetName()).c_str(), "L");
		}
	}
	/*
	if(BkgnSpect){
		if(nDs>0){
			BkgnSpect->Draw("same");
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}else{
			BkgnSpect->Draw();
			leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
		}
	}
	*/
	
	if(BkgnSpect){
		leg->AddEntry(BkgnSpect->GetHisto(), (BkgnSpect->GetName()).c_str(), "L");
	}
	
	leg->Draw();
}


#include "../source/src/loadSPE.cc"
#include "../source/src/convert_histo_ENR.cc"