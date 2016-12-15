#include <vector>
#include <string>
#include <sstream>

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
#include "misc.h"
#include "GatorCalibFunc.h"
//#include "GatorLoadData.h"


using namespace std;


int GatorCalibScriptBAT(string calibset, string configfile, bool recreate)
{
	const string GatorSystem = getEnvVar("GATOR_SYS");
	if(GatorSystem==string("")){
		cerr << "\n\nERROR --> Environment variable \"GATOR_SYS\" is not set cannot continue.\n" << endl;
		return -1;
	}
	
	//Search for the file containing the path to the calibration archive
	string archivepathfile;
	if( GatorSystem.at(GatorSystem.length()-1) != '/' ){
		archivepathfile = GatorSystem + string("/etc/CalibArchive.conf");
	}else{
		archivepathfile = GatorSystem + string("etc/CalibArchive.conf");
	}
	
	string archivedir("");
	if(!fexist(archivepathfile)){
		cerr << "\n\nERROR --> Cannot open/find file <" << archivepathfile << ">\n" << endl;
		return -2;
	}else{
		stringstream line; line.str("");
		ifstream infile(archivepathfile.c_str());
		if(getline(infile, archivedir)){
			line << archivedir;
			line >> archivedir;
		}

		if(archivedir.at(archivedir.length()-1) != '/' ) archivedir = archivedir + string("/");
	}
	
	
	string calibdir = archivedir + calibset + string("/");
	
	//Load the spectrum of the calibs
	double calibtime = 0.;
	cout << "Calibration dir: <" << calibdir << ">" << endl;
	TH1D* MCAhisto = loadSpe(calibdir.c_str(),calibtime);
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	//Initialize the lines for the calibration from the "configfile"
	cout << "Configuration file: <" << configfile << ">" << endl;
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
	
	TTree *t1=NULL, *t1_old=NULL;
	if(!recreate){
		t1_old = (TTree*)outfile->Get("linestree");
		if(!t1_old){
			recreate = true;//Behaves like if the tree was recreated even if the file is opened in update mode
		}
	}
	
	t1 = new TTree("linestree","Results of the fits for each line");
	
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
	
	TH1D* p_line_histo = NULL;
	
	
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
	
	t1 -> Branch("histo", "TH1D", &p_line_histo);
	
	
	if(!recreate){
		t1_old -> SetBranchAddress("massNum",&massN);
		t1_old -> SetBranchAddress("element",&pelement);
		t1_old -> SetBranchAddress("litEn",&litEn);
		t1_old -> SetBranchAddress("litEn_err",&litEn_err);
		t1_old -> SetBranchAddress("mean",&mean);
		t1_old -> SetBranchAddress("mean_err",&mean_err);
		t1_old -> SetBranchAddress("sigma",&sigma);
		t1_old -> SetBranchAddress("sigma_err",&sigma_err);
		t1_old -> SetBranchAddress("beta",&beta);
		t1_old -> SetBranchAddress("beta_err",&beta_err);
		t1_old -> SetBranchAddress("ampl",&ampl);
		t1_old -> SetBranchAddress("ampl_err",&ampl_err);
		t1_old -> SetBranchAddress("tail",&tail);
		t1_old -> SetBranchAddress("tail_err",&tail_err);
		//t1_old -> SetBranchAddress("ratio",&ratio);
		//t1_old -> SetBranchAddress("ratio_err",&ratio_err);
		t1_old -> SetBranchAddress("step",&step);
		t1_old -> SetBranchAddress("step_err",&step_err);
		t1_old -> SetBranchAddress("cost",&cost);
		t1_old -> SetBranchAddress("cost_err",&cost_err);
		t1_old -> SetBranchAddress("p_value",&p_value);
		t1_old -> SetBranchAddress("p_value_ndof",&p_value_ndof);
		
		t1_old -> SetBranchAddress("histo", &p_line_histo);
	}
	
	
	
	for(int iLine=0; iLine<nlines; iLine++){
		string ans("");
		while(true){
			cout << "\nDo you want to fit the " << CalibLinesVec[iLine].massN << "-" << CalibLinesVec[iLine].element << " line at " << CalibLinesVec[iLine].litEn << " keV line? [y/n]\n" << "> ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("N") || ans!=string("n") || ans!=string("N")) break;
		}
		
		if(ans==string("y")  || ans==string("Y")){
			if(doFitBAT(MCAhisto,CalibLinesVec[iLine])){
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
				
				p_line_histo = CalibLinesVec[iLine].histo;
				
				t1->Fill();
				t1->AutoSave();
			}else{
				if(t1_old){
					for(int iEnt=0; iEnt<t1_old->GetEntries(); iEnt++){
						t1_old->GetEntry(iEnt);
						if( (massN==CalibLinesVec[iLine].massN) && (element==CalibLinesVec[iLine].element) && (litEn==CalibLinesVec[iLine].litEn)){
							t1->Fill();
							t1->AutoSave();
							break;
						}
					}
				}
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



int GatorCalibScriptLL(string calibset, string configfile, bool recreate)
{
	const string GatorSystem = getEnvVar("GATOR_SYS");
	if(GatorSystem==string("")){
		cerr << "\n\nERROR --> Environment variable \"GATOR_SYS\" is not set cannot continue.\n" << endl;
		return -1;
	}
	
	//Search for the file containing the path to the calibration archive
	string archivepathfile;
	if( GatorSystem.at(GatorSystem.length()-1) != '/' ){
		archivepathfile = GatorSystem + string("/etc/CalibArchive.conf");
	}else{
		archivepathfile = GatorSystem + string("etc/CalibArchive.conf");
	}
	
	string archivedir("");
	if(!fexist(archivepathfile)){
		cerr << "\n\nERROR --> Cannot open/find file <" << archivepathfile << ">\n" << endl;
		return -2;
	}else{
		stringstream line; line.str("");
		ifstream infile(archivepathfile.c_str());
		if(getline(infile, archivedir)){
			line << archivedir;
			line >> archivedir;
		}

		if(archivedir.at(archivedir.length()-1) != '/' ) archivedir = archivedir + string("/");
	}
	
	
	string calibdir = archivedir + calibset + string("/");
	
	//Load the spectrum of the calibs
	double calibtime = 0.;
	cout << "Calibration dir: <" << calibdir << ">" << endl;
	TH1D* MCAhisto = loadSpe(calibdir.c_str(),calibtime);
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	//Initialize the lines for the calibration from the "configfile"
	cout << "Configuration file: <" << configfile << ">" << endl;
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
	
	TTree *t1=NULL, *t1_old=NULL;
	if(!recreate){
		t1_old = (TTree*)outfile->Get("linestree");
		if(!t1_old){
			recreate = true;//Behaves like if the tree was recreated even if the file is opened in update mode
		}
	}
	
	t1 = new TTree("linestree","Results of the fits for each line");
	
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
	
	double chi2;
	double chi2ndof;
	
	double p_value;
	double p_value_ndof;
	
	TH1D* p_line_histo = NULL;
	
	
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
	
	t1 -> Branch("chi2",&chi2,"chi2/D");
	t1 -> Branch("chi2ndof",&chi2ndof,"chi2ndof/D");
	
	t1 -> Branch("p_value",&p_value,"p_value/D");
	t1 -> Branch("p_value_ndof",&p_value_ndof,"p_value_ndof/D");
	
	t1 -> Branch("histo", "TH1D", &p_line_histo);
	
	
	if(!recreate){
		t1_old -> SetBranchAddress("massNum",&massN);
		t1_old -> SetBranchAddress("element",&pelement);
		t1_old -> SetBranchAddress("litEn",&litEn);
		t1_old -> SetBranchAddress("litEn_err",&litEn_err);
		t1_old -> SetBranchAddress("mean",&mean);
		t1_old -> SetBranchAddress("mean_err",&mean_err);
		t1_old -> SetBranchAddress("sigma",&sigma);
		t1_old -> SetBranchAddress("sigma_err",&sigma_err);
		t1_old -> SetBranchAddress("beta",&beta);
		t1_old -> SetBranchAddress("beta_err",&beta_err);
		t1_old -> SetBranchAddress("ampl",&ampl);
		t1_old -> SetBranchAddress("ampl_err",&ampl_err);
		t1_old -> SetBranchAddress("tail",&tail);
		t1_old -> SetBranchAddress("tail_err",&tail_err);
		//t1_old -> SetBranchAddress("ratio",&ratio);
		//t1_old -> SetBranchAddress("ratio_err",&ratio_err);
		t1_old -> SetBranchAddress("step",&step);
		t1_old -> SetBranchAddress("step_err",&step_err);
		t1_old -> SetBranchAddress("cost",&cost);
		t1_old -> SetBranchAddress("cost_err",&cost_err);
		
		if( t1_old->GetBranch("chi2") ) t1_old -> SetBranchAddress("chi2", &chi2);
		if( t1_old->GetBranch("chi2ndof") ) t1_old -> SetBranchAddress("chi2ndof", &chi2ndof);
		
		t1_old -> SetBranchAddress("p_value",&p_value);
		t1_old -> SetBranchAddress("p_value_ndof",&p_value_ndof);
		if( t1_old->GetBranch("histo") ) t1_old -> SetBranchAddress("histo", &p_line_histo);
	}
	
	
	
	for(int iLine=0; iLine<nlines; iLine++){
		string ans("");
		while(true){
			cout << "\nDo you want to fit the " << CalibLinesVec[iLine].massN << "-" << CalibLinesVec[iLine].element << " line at " << CalibLinesVec[iLine].litEn << " keV line? [y/n]\n" << "> ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("N") || ans!=string("n") || ans!=string("N")) break;
		}
		
		if(ans==string("y")  || ans==string("Y")){
			if(doFitLL(MCAhisto,CalibLinesVec[iLine])){
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
				step = CalibLinesVec[iLine].step;
				step_err = CalibLinesVec[iLine].step_err;
				cost = CalibLinesVec[iLine].cost;
				cost_err = CalibLinesVec[iLine].cost_err;
				
				chi2 = CalibLinesVec[iLine].chi2;
				chi2ndof = CalibLinesVec[iLine].chi2ndof;
				
				p_value = CalibLinesVec[iLine].p_value;
				p_value_ndof = CalibLinesVec[iLine].p_value_ndof;
				
				p_line_histo = CalibLinesVec[iLine].histo;
				
				t1->Fill();
				t1->AutoSave();
			}else{
				if(t1_old){
					for(int iEnt=0; iEnt<t1_old->GetEntries(); iEnt++){
						t1_old->GetEntry(iEnt);
						if( (massN==CalibLinesVec[iLine].massN) && (element==CalibLinesVec[iLine].element) && (litEn==CalibLinesVec[iLine].litEn)){
							t1->Fill();
							t1->AutoSave();
							break;
						}
					}
				}
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
