#include <vector>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#include <BAT/BCAux.h>
//#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>

#include "screenfncs.h"
#include "GatorStructs.h"
#include "GatorCalibFunc.h"
//#include "GatorLoadData.h"


using namespace std;

int GatorCalibScript(string calibset, string configfile, bool recreate){
	
	//const string archivedir("/home/atp/fpiastra/Gator/calibrations/archive/");
	const string archivedir("/Users/francesco/PhD/Gator/data_analysis/calibrations/archive/");
	
	string calibdir = archivedir + calibset + string("/");
	
	//Load the spectrum of the calibs
	double calibtime = 0.;
	TH1F* MCAhisto = loadSpe(calibdir.c_str(),calibtime);
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	//Initialize the lines for the calibration from the "configfile"
	vector<string> linesinfo = loadConfFile(configfile);
	vector<CalibLine> CalibLinesVec;
	int nlines = linesinfo.size();
	for(int iLine=0; iLine<nlines; iLine++){
		CalibLinesVec.push_back(lineInit(linesinfo[iLine]));
	}
	
	//Generate a root file and a TTree with the fitted lines values
	string openmode;
	if(recreate){
		openmode = string("recreate");
	}else{
		openmode = string("update");//Default mode
	}
	
	
	gROOT ->ProcessLine("#include <string>");
	
	TFile *outfile = new TFile((calibdir+string("fittedlines.root")).c_str(),openmode.c_str());
	outfile->cd();
	
	bool newfile = false;
	
	TTree *t1;
	if(recreate){
		t1 = new TTree("linestree","Results of the fits for each line");
		t1->SetDirectory(outfile);
	}else{
		t1 = (TTree*)outfile->Get("linestree");
		if(!t1){
			t1 = new TTree("linestree","Results of the fits for each line");
			t1->SetDirectory(outfile);
			newfile = true;
		}
	}
	
	//Variables to put in tree
	string massN;
	string *pmassN = &massN;
	string element;
	string *pelement = &element;
	double litEn;
	double litEn_err;
	double mean;
	double mean_err;
	double sigma;
	double sigma_err;
	double beta;
	double beta_err;
	double ampl;
	double ampl_err;
	double tail;
	double tail_err;
	//double ratio;
	//double ratio_err;
	double step;
	double step_err;
	double cost;
	double cost_err;
	double p_value;
	double p_value_ndof;
	
	if(recreate || newfile){
		t1 -> Branch("massNum","string",&massN);
		t1 -> Branch("element","string",&pelement);
		t1 -> Branch("litEn",&litEn,"litEn/D");
		t1 -> Branch("litEn_err",&litEn_err,"litEn_err/D");
		t1 -> Branch("mean",&mean,"mean/D");
		t1 -> Branch("mean_err",&mean_err,"mean_err/D");
		t1 -> Branch("sigma",&sigma,"sigma/D");
		t1 -> Branch("sigma_err",&sigma_err,"sigma_err/D");
		t1 -> Branch("beta",&beta,"beta/D");
		t1 -> Branch("beta_err",&beta_err,"beta_err/D");
		t1 -> Branch("ampl",&ampl,"ampl/D");
		t1 -> Branch("ampl_err",&ampl_err,"ampl_err/D");
		t1 -> Branch("tail",&tail,"tail/D");
		t1 -> Branch("tail_err",&tail_err,"tail_err/D");
		//t1 -> Branch("ratio",&ratio,"ratio/D");
		//t1 -> Branch("ratio_err",&ratio_err,"ratio_err/D");
		t1 -> Branch("step",&step,"step/D");
		t1 -> Branch("step_err",&step_err,"step_err/D");
		t1 -> Branch("cost",&cost,"cost/D");
		t1 -> Branch("cost_err",&cost_err,"cost_err/D");
		t1 -> Branch("p_value",&p_value,"p_value/D");
		t1 -> Branch("p_value_ndof",&p_value_ndof,"p_value_ndof/D");
	}else{
		t1 -> SetBranchAddress("massNum",&massN);
		t1 -> SetBranchAddress("element",&pelement);
		t1 -> SetBranchAddress("litEn",&litEn);
		t1 -> SetBranchAddress("litEn_err",&litEn_err);
		t1 -> SetBranchAddress("mean",&mean);
		t1 -> SetBranchAddress("mean_err",&mean_err);
		t1 -> SetBranchAddress("sigma",&sigma);
		t1 -> SetBranchAddress("sigma_err",&sigma_err);
		t1 -> SetBranchAddress("beta",&beta);
		t1 -> SetBranchAddress("beta_err",&beta_err);
		t1 -> SetBranchAddress("ampl",&ampl);
		t1 -> SetBranchAddress("ampl_err",&ampl_err);
		t1 -> SetBranchAddress("tail",&tail);
		t1 -> SetBranchAddress("tail_err",&tail_err);
		//t1 -> SetBranchAddress("ratio",&ratio);
		//t1 -> SetBranchAddress("ratio_err",&ratio_err);
		t1 -> SetBranchAddress("step",&step);
		t1 -> SetBranchAddress("step_err",&step_err);
		t1 -> SetBranchAddress("cost",&cost);
		t1 -> SetBranchAddress("cost_err",&cost_err);
		t1 -> SetBranchAddress("p_value",&p_value);
		t1 -> SetBranchAddress("p_value_ndof",&p_value_ndof);
	}
	
	
	
	
	for(int iLine=0; iLine<nlines; iLine++){
		string ans("");
		while(true){
			cout << "\nDo you want to fit the " << CalibLinesVec[iLine].massN << "-" << CalibLinesVec[iLine].element << " line at " << CalibLinesVec[iLine].litEn << " keV line? [y/n]\n" << "> ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("N") || ans!=string("n") || ans!=string("N")) break;
		}
		
		if(ans==string("y")  || ans==string("Y")){
			if(doFit(MCAhisto,CalibLinesVec[iLine])){
				massN = CalibLinesVec[iLine].massN;
				element = CalibLinesVec[iLine].element;
				litEn = CalibLinesVec[iLine].litEn;
				litEn_err = CalibLinesVec[iLine].litEn_err;
				mean = CalibLinesVec[iLine].mean;
				mean_err = CalibLinesVec[iLine].mean_err;
				sigma = CalibLinesVec[iLine].sigma;
				sigma_err = CalibLinesVec[iLine].sigma_err;
				beta = CalibLinesVec[iLine].beta;
				beta_err = CalibLinesVec[iLine].beta_err;
				ampl = CalibLinesVec[iLine].ampl;
				ampl_err = CalibLinesVec[iLine].ampl_err;
				tail = CalibLinesVec[iLine].tail;
				tail_err = CalibLinesVec[iLine].tail_err;
				//ratio = CalibLinesVec[iLine].ratio;
				//ratio_err = CalibLinesVec[iLine].ratio_err;
				step = CalibLinesVec[iLine].step;
				step_err = CalibLinesVec[iLine].step_err;
				cost = CalibLinesVec[iLine].cost;
				cost_err = CalibLinesVec[iLine].cost_err;
				p_value = CalibLinesVec[iLine].p_value;
				p_value_ndof = CalibLinesVec[iLine].p_value_ndof;
				
				t1->Fill();
				t1->AutoSave();
			}
		}
	}
	
	if(t1) t1->Write();
	if(outfile){
		outfile->Close();
		delete outfile;
	}
	
	return 0;
	
}
