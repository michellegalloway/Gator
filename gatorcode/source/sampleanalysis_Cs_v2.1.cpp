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

//////////////////////////////////////////////////////////////////////////////////////////

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


void sampleanalysis(){
	
	gStyle->SetOptStat("");
	c1= new TCanvas("c1","");
	c1->cd();
	
	//The following lines are for closing the application when the canvas c1 is closed.

#ifndef __CINT__
	c1 -> Connect("TCanvas", "Closed()", "TApplication", gApplication,"Terminate()");
#endif

	
//////////////////////////////////////////////////////////////////////////////////////////		
// hard-coded input variables:
	
	string datadir = "/home/atp/galloway/Gator_Screening/Data/Samples/PMT_R11410/PMT_15_20150404/SPE/"; // This could be BG SPE dir if "BGswitch" is false	
	string workdir = "/home/atp/galloway/Gator_Screening/Data/Samples/PMT_R11410/PMT_15_20150404/"; // output files directory
	

	//Cs line (only) reevaluated following contamination, then recalibration
	string backgroundSPEdir = "/home/atp/galloway/Gator_Screening/Data/Background_PMTs/PTFE_holders_Oct2014/SPE/"; // For PMTs this is the measurements with the PTFE holders	
	string calibdir = "/home/atp/galloway/Gator_Screening/Calibrations/bkgrd_Cs_201410/"; // detector recalibrated Nov 1, 2014
	string calibdir_bg = "/home/atp/galloway/Gator_Screening/Calibrations/bkgrd_all_201407/";	
	string configdir = "/home/atp/galloway/Gator_Screening/Data/Samples/PMT_R11410/PMT_15_20141210/";
	
string linesfilename = configdir + string("137Cslines.list"); //It doesn't matter if the "Eff_rootfile" switch is true. In such a case it is not used at all.
	

	string calfilename = "calibration.root";
	calfilename = calibdir + calfilename;
	

	string calfilename_bg = "calibration.root";
	calfilename_bg = calibdir_bg + calfilename_bg;
	
	
	string outfilename = "PMTs_Cs.txt";
	outfilename = workdir + outfilename;
	
	string efftabfilename = "PMTs_eff_Cs.txt";
	efftabfilename = workdir + efftabfilename;
	
	string savefilename = "PMTs_Cs.root";
	savefilename = workdir + savefilename;
	
	Double_t quantity = 15; //This is the sample weight in kg!!! Or can even be the number of the pieces (PMTs e.g.)

// HARD-CODED input stuff is finished.

//////////////////////////////////////////////////////////////////////////////////////////	


	
	GatorSampleAnalysis analyser(workdir);
	
	
	
	analyser.SetSampleCalib(calibdir);
	analyser.SetBackgroundCalib(calibdir_bg);
	
	if(!analyser.LoadData(datadir)){
		exit(-1);
	}
	
	if(!analyser.LoadBackground(backgroundSPEdir)){
		exit(-1);
	}
	if(!analyser.LoadLines(linesfilename)){
		exit(-1);
	}
	//return;
	
	
	//analyser.SetLeftCompRegion(651, 652.8, 0);
	//analyser.SetRightCompRegion(663, 664.8, 0);
	
	analyser.DefaultIntervals();
	
	analyser.SetSampleQuntity(quantity);
	analyser.SetQuantityUnit("PMT");
	
	
	
	//Start the analysis here
	
	analyser.WriteEffTable(efftabfilename);
	analyser.WriteOutputTable(outfilename);
	
	analyser.DrawActivityPlots(c1);

	
	theApp -> Run(kTRUE);
	
	
	return;
}
//////////////////////////////////////////////////////////////////
//--------------------------------------------------------------//
//                SCRIPT FUNCTION FINISHED                      //
//--------------------------------------------------------------//
//////////////////////////////////////////////////////////////////
