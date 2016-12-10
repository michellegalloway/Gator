#include "GatorTreeDataLoader.hh"
#include "GatorParam.hh"

#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"


using namespace std;

GatorTreeDataLoader::GatorTreeDataLoader():
fSourceFileName(""),
fSigLevel(0),
fTreeDataLoaderActive(false),
fDataLoaded(false)
{
}


GatorTreeDataLoader::GatorTreeDataLoader( string sourcefilename, string treename ){
	
	fSourceFileName = sourcefilename;
	fTreeName = treename;
	fTreeDataLoaderActive = true;
	fSize = 0;
	fSigLevel = 0;
	fTreeEntries = 0;
	fDataLoaded = false;
	fSingleLine = false;
	
	
	gROOT->ProcessLine("#include <vector>");
	
	int errcode = LoadTree();
	if(errcode>0){
		cerr << "\n\nGatorTreeDataLoader --> ERROR: Problem in reding TParameters in <" << fSourceFileName << ">.\nAborting program with error code " << errcode << "\n" << endl;
		exit(errcode);
	}
	if(errcode<0){
		cerr << "\n\nGatorTreeDataLoader --> ERROR: Generic problem in <" << fSourceFileName << ">.\nAborting program with error code " << errcode << "\n" << endl;
		exit(errcode);
	}
	
}


int GatorTreeDataLoader::LoadTree(){
	
	cout << "\nEntering in GatorTreeDataLoader::LoadTree()\n" << endl;
	
	stringstream ss_tmp; //Always useful
	
	//Determine how many lines are in the dataset
	
	TFile *pSourceFile = new TFile(fSourceFileName.c_str(),"read");
	
	if(!pSourceFile) return(1);
	
	TParameter<int> *tp_nLines = (TParameter<int>*)pSourceFile->Get("nLines");
	
	if(!tp_nLines){
		pSourceFile -> Close();
		return(2);
	}
	m_nLines = (unsigned)tp_nLines->GetVal();
	cout << "Number of lines: " << m_nLines << endl;
	
	if(m_nLines>0){
		
		if(m_nLines==1){
			fSingleLine = true;
		}
		
		//m_nPar = 2*m_nLines +1;
		
	}else{
		pSourceFile -> Close();
		return(3);
	}
	
	//Determine the lenght of the spectra of the lines (for all the lines is the same)
	TParameter<int> *tp_sizespec = (TParameter<int>*)pSourceFile->Get("sizespec");
	
	if(!tp_sizespec){
		pSourceFile -> Close();
		return(2);
	}
	fSize = (unsigned)tp_sizespec->GetVal();
	cout << "Each line's spectrum is long " << fSize << " MCA channels" << endl;
	
	//Determine the signal from tParameter
	TParameter<double> *tp_sigLevel = (TParameter<double>*)pSourceFile->Get("siglevel");
	
	if(!tp_sigLevel){
		pSourceFile -> Close();
		return(7);
	}
	fSigLevel = tp_sigLevel->GetVal();
	cout << "Nominal level of the signal: " << fSigLevel << " cnts" << endl;
	
	TParameter<double> *tp_snrNomin = (TParameter<double >*)pSourceFile->Get("nomsnr");
	if(!tp_snrNomin){
		pSourceFile -> Close();
		return(8);
	}
	fSnrLevel = tp_snrNomin->GetVal();
	cout << "Nominal SNR value (only on first line): " << fSnrLevel << endl;
	
	
	//Determine the sigma from tParameter
	vector<double> *v_sigma = (vector<double> *)pSourceFile->Get("sigma");
	
	if(!v_sigma){
		pSourceFile -> Close();
		return(4);
	}
	fSigma = *(v_sigma);
	cout << "Sigma value of the first line " << fSigma.at(0) << " MCA channels" << endl;
	
	
	//Determine the background level from tParameter
	vector<double> *v_bglevel = (vector<double >*)pSourceFile->Get("bglevel");
	
	if(!v_bglevel){
		pSourceFile -> Close();
		return(5);
	}
	fBgLevel = *(v_bglevel);
	cout << "Background level for the first line: " << fBgLevel.at(0) << " cnts/ch" << endl;
	
	//Determine the efficiencies for each line
	vector<double> *v_effic = (vector<double >*)pSourceFile->Get("effic");
	
	if(!v_effic){
		pSourceFile -> Close();
		return(6);
	}
	fEffic = *(v_effic);
	cout << "Efficiency for the first line: " << fEffic.at(0) << endl;
	
	
	TTree *t1 = (TTree*)pSourceFile->Get(fTreeName.c_str());
	cout << "Searching for the tree \"" << fTreeName << "\"...  ";
	if(!t1){
		pSourceFile -> Close();
		return(-1);
	}
	fTreeEntries = t1->GetEntries();
	
	cout << "found " << fTreeEntries << " entries." << endl;
	
	pSourceFile -> Close();
	
	fDataLoaded = true;
	
	cout << "\nExiting from GatorTreeDataLoader::LoadTree() function" << endl << endl;
	
	return(0);
}
