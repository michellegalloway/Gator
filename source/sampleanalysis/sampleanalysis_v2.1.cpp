#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
//#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>
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
#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
#include <TFitResult.h>

#include "GatorStructs.h"
#include "GatorCounter.hh"
#include "GatorDataLoader.hh"
#include "GatorSampleAnalysis.hh"
#include "screenfncs.h"


using namespace std;

void sampleanalysis();
void wait();

TCanvas* c1;

extern TApplication* theApp=NULL;

TApplication *myApp;

int main(){
	
	if(!theApp) theApp = new TApplication("theApp",0,0);
	sampleanalysis();
	
	return 0;
}


void wait()
{
	c1 -> Update(); c1 -> WaitPrimitive();
}
//////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
////////////---------- SCRIPT FUNCTION ------------//////////////
/////////////////////////////////////////////////////////////////


void sampleanalysis(){
	
	gStyle->SetOptStat("");
	
	c1= new TCanvas("c1","");
	c1->cd();
	
	
	//The following lines are for closing the application when the canvas c1 is closed.

#ifndef __CINT__
	c1 -> Connect("TCanvas", "Closed()", "TApplication", gApplication, "Terminate()");
#endif

	
	string GatorSystem = getEnvVar("GATOR_SYS");
	if(GatorSystem==string("")){
		cerr << "\n\nERROR --> Environment variable \"GATOR_SYS\" is not set cannot continue.\n" << endl;
		return;
	}
	string archivepathfile;
	if( GatorSystem.at(GatorSystem.length()-1) != '/' ) GatorSystem = GatorSystem + string("/");
	
	
// HARD-CODED input variables:
	
	string workdir = GatorSystem + "screening/SAMPLE_NAME/"; //Where the output files are placed
	
	//Directory where are the SPE files relative to the sample 
	string datadir = GatorSystem + "screening/SAMPLE_NAME/SPE/"; //This could be BG SPE dir if "BGswitch" is false
	
	//Directory where are the SPE files relative to the background 
	string backgroundSPEdir = GatorSystem + "background/archive/2014/";
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/analysis/samples/PTFE_holders/SPE2/"; //This could be the directory of screening files done before activation (the Xe activation work for example)
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/analysis/sessions/BG_137Cs/EmptyCavity/"; //6days. Keep here otherwise I forget!!!
	//string backgroundSPEdir = "/Users/francesco/PhD/Gator/analysis/sessions/BG_137Cs/Cs+PMTholders/"; //6.8days. Keep here otherwise I forget!!!
	
	string configdir = workdir;
	string linesfilename = configdir + string("lines.list");
	
	string calibdir = GatorSystem + "calibrations/archive/2013.07.01/"; //This is relative to the sample
	string calfilename = "calibration.root"; //Don't touch
	calfilename = calibdir + calfilename;
	
	string calibdir_bg = GatorSystem + "calibrations/archive/2013.12.23/"; //This is relative to the background
	//string calibdir_bg = "/Users/francesco/PhD/Gator/calibrations/archive/2014.10.10/";
	string calfilename_bg = "calibration.root"; //Don't touch
	calfilename_bg = calibdir_bg + calfilename_bg;
	
	
	//Output files:
	
	//Results table in wiki format
	string outfilename = "CeramicStems_R11410-21_STD.txt"; //Edit this
	outfilename = workdir + outfilename;
	
	//Table of BRs and efficiencies in wiki format
	string efftabfilename = "CeramicStems_R11410-21_STD_efftab.txt"; //Edit this
	efftabfilename = workdir + efftabfilename;
	
	//Everithing in the root format
	string savefilename = "CeramicStems_R11410-21_STD.root";
	savefilename = workdir + savefilename;
	
	
	Double_t quantity = 6; //This is the sample weight in kg!!! Or can even be the number of the pieces (PMTs e.g.)

// HARD-CODED input stuff is finished.
	
	
	
	
	
	GatorSampleAnalysis *analyser = new GatorSampleAnalysis(workdir);
	
	
	
	analyser->SetSampleCalib(calibdir);
	analyser->SetBackgroundCalib(calibdir_bg);
	
	if(!(analyser->LoadData(datadir))){
		cerr << "\nERROR in loading data. EXIT!" << endl;
		exit(-1);
	}
	
	if(!(analyser->LoadBackground(backgroundSPEdir))){
		cerr << "\nERROR in loading background. EXIT!" << endl;
		exit(-1);
	}
	if(!(analyser->LoadLines(linesfilename))){
		cerr << "\nERROR in loading the list of lines. EXIT!" << endl;
		exit(-1);
	}
	//return;
	
	//analyser->ForcePeakBg();
	//analyser->SetLeftCompRegion(654.076, 655.918, 0);
	//analyser->SetRightCompRegion(663.502, 665.347, 0);
	
	analyser->DefaultIntervals();
	
	analyser->SetSampleQuntity(quantity);
	analyser->SetQuantityUnit("PMT"); //Edit this
	
	
	
	//Start the analysis here
	
	analyser->WriteEffTable(efftabfilename);
	analyser->WriteOutputTable(outfilename);
	
	analyser->DrawActivityPlots(c1); //Not working at the moment.

	
	theApp -> Run(kTRUE);
	
	
	return;
}
//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//                SCRIPT FUNCTION FINISHED                      //
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
