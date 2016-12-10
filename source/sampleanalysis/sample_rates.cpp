#include <cstdlib>
#include <unistd.h>

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

#include <TROOT.h>
//#include <TSystem.h>
//#include <TStyle.h>
//#include <TFile.h>
//#include <TTree.h>
#include <TH1F.h>
//#include <TH2F.h>
//#include <TH1D.h>
//#include <TH2D.h>
//#include <TF1.h>
//#include <TLegend.h>
//#include <TCut.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
//#include <TIterator.h>
//#include <TList.h>
//#include <TMultiGraph.h>
#include <TMath.h>
#include <TApplication.h>
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
#include "trigrate.hh"

using namespace std;


int main(int argc, char** argv){
	
	string SPEdir;
	string OutDir;
	string OutFile;
	stringstream StartEnrStr("");
	stringstream EndEnrStr("");
	Double_t StartCh, EndCh;
	
	//Read the options and set the flags and filenames
	//string sourceName;
	//string outFileName;
	bool flag_setSPEDir = false;
	bool flag_setOutDir = false;
	bool flag_setOutFile = false;
	bool flag_setStartCh = false;
	bool flag_setEndCh = false;
	char c=0;
	while((c=getopt(argc,argv,"i:o:s:e:f:"))!=-1){
		switch(c)
		{
			case 'i':
				SPEdir = optarg;
				flag_setSPEDir = true;
				break;
			case 'o':
				OutDir = optarg;
				flag_setOutDir = true;
				break;
			case 's':
				StartEnrStr << optarg;
				if(!(StartEnrStr>>StartCh)){
					cout << endl << "ERROR:: You should plug in a numeric value for the start channel" << endl;
					return -1;
				}else{
					flag_setStartCh = true;
				}
				break;
			case 'e':
				EndEnrStr << optarg;
				if(!(EndEnrStr>>EndCh)){
					cout << endl << "ERROR:: You should plug in a numeric value for the end channel" << endl;
					return -1;
				}else{
					flag_setEndCh = true;
				}
				break;
			case 'f':
				OutFile = optarg;
				flag_setOutFile = true;
				break;
			default:
				cout << "Unrecognized option: " << c << endl;
				return -1;
		};
	};
	
	if(!flag_setSPEDir){
		cout << "ERROR:: SPE directory was not set." << endl;
		cout << endl << "Usage: sample_rate -i <SPEdir> [-o <outputdir>|-h <ADC start counting>]" << endl << endl;
	}else{
		if(SPEdir.substr(SPEdir.size()-1,1)!=string("/")) SPEdir = SPEdir + string("/");
	}
	
	if(!flag_setOutDir){
		OutDir = SPEdir + string("../");
	}else{
		if(OutDir.substr(OutDir.size()-1,1)!=string("/")) OutDir = OutDir + string("/");
	}
	
	cout << "SPE directory: <" << SPEdir << ">\n\n" << endl;
	
	trigrate* pTrigRate = new trigrate();
	
	if(flag_setStartCh&&flag_setEndCh){
		pTrigRate->FindSPE(SPEdir.c_str(), StartCh, EndCh, true);
	}else{
		pTrigRate->FindSPE(SPEdir.c_str(),true);
	}
	
	TGraphErrors* gr_rates = pTrigRate->MakeFullPlot();
	
	TCanvas* c1 = new TCanvas("c1");
	
	gr_rates->Draw("ap*");
	
	if(flag_setOutFile){
		OutFile = OutFile + string(".root");
		c1->SaveAs((OutDir+OutFile).c_str());
	}else{
		c1->SaveAs((OutDir+string("trigrate_plot.root")).c_str());
	}
	
	return 0;		
}