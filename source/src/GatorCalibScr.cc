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

int GatorCalibScriptBAT(string calibset, string configfile, bool recreate, bool update)
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
	
	if(recreate && update){
		cerr << "\n\nERROR --> Both the \"recreate\" and \"update\" modes are set! Cannot use conficting options!\n" << endl;
		return -3;
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
		openmode = string("update"); //Default mode
	}
	
	
	gROOT ->ProcessLine("#include <string>");
	
	TFile *outfile = new TFile((calibdir+string("fittedlines.root")).c_str(),openmode.c_str());
	outfile->cd();
	
	bool newfile = false;
	
	
	TTree *t1=NULL, *t1_old=NULL;
	if(!recreate){
		t1_old = (TTree*)outfile->Get("linestree");
		if(!t1_old){
			recreate = true; //Behaves like if the tree was recreated even if the file is opened in update mode
		}
	}
	
	t1 = new TTree("linestree","Results of the fits for each line");
	
	
	CalibLine LinkedLineStruct; //This is to address the trees.
	
	//Used only to fulfil the branching address syntax
	string *p_massN = &LinkedLineStruct.massN;
	string *p_element = &LinkedLineStruct.element;
	
	t1 -> Branch("massNum","string",&p_massN);
	t1 -> Branch("element","string",&p_element);
	t1 -> Branch("litEn",&LinkedLineStruct.litEn,"litEn/D");
	t1 -> Branch("litEn_err",&LinkedLineStruct.litEn_err,"litEn_err/D");
	t1 -> Branch("mean",&LinkedLineStruct.mean,"mean/D");
	t1 -> Branch("mean_err",&LinkedLineStruct.mean_err,"mean_err/D");
	t1 -> Branch("sigma",&LinkedLineStruct.sigma,"sigma/D");
	t1 -> Branch("sigma_err",&LinkedLineStruct.sigma_err,"sigma_err/D");
	t1 -> Branch("beta",&LinkedLineStruct.beta,"beta/D");
	t1 -> Branch("beta_err",&LinkedLineStruct.beta_err,"beta_err/D");
	t1 -> Branch("ampl",&LinkedLineStruct.ampl,"ampl/D");
	t1 -> Branch("ampl_err",&LinkedLineStruct.ampl_err,"ampl_err/D");
	t1 -> Branch("tail",&LinkedLineStruct.tail,"tail/D");
	t1 -> Branch("tail_err",&LinkedLineStruct.tail_err,"tail_err/D");
	//t1 -> Branch("ratio",&LinkedLineStruct.ratio,"ratio/D");
	//t1 -> Branch("ratio_err",&LinkedLineStruct.ratio_err,"ratio_err/D");
	t1 -> Branch("step",&LinkedLineStruct.step,"step/D");
	t1 -> Branch("step_err",&LinkedLineStruct.step_err,"step_err/D");
	t1 -> Branch("cost",&LinkedLineStruct.cost,"cost/D");
	t1 -> Branch("cost_err",&LinkedLineStruct.cost_err,"cost_err/D");
	
	t1 -> Branch("chi2",&LinkedLineStruct.chi2,"chi2/D");
	t1 -> Branch("chi2ndof",&LinkedLineStruct.chi2ndof,"chi2ndof/D");
	
	t1 -> Branch("p_value",&LinkedLineStruct.p_value,"p_value/D");
	t1 -> Branch("p_value_ndof",&LinkedLineStruct.p_value_ndof,"p_value_ndof/D");
	
	t1 -> Branch("histo", "TH1D", &LinkedLineStruct.histo, 32000, 0);
	t1 -> Branch("fit", "TF1", &LinkedLineStruct.fit, 32000, 0);
	
	
	vector<CalibLine> OldCalibLinesVector; //Used only in the forced update mode
	
	if(!recreate){
		if(t1_old){
			
			
			t1_old -> SetBranchAddress("massNum",&p_massN);
			t1_old -> SetBranchAddress("element",&p_element);
			t1_old -> SetBranchAddress("litEn",&LinkedLineStruct.litEn);
			t1_old -> SetBranchAddress("litEn_err",&LinkedLineStruct.litEn_err);
			t1_old -> SetBranchAddress("mean",&LinkedLineStruct.mean);
			t1_old -> SetBranchAddress("mean_err",&LinkedLineStruct.mean_err);
			t1_old -> SetBranchAddress("sigma",&LinkedLineStruct.sigma);
			t1_old -> SetBranchAddress("sigma_err",&LinkedLineStruct.sigma_err);
			t1_old -> SetBranchAddress("beta",&LinkedLineStruct.beta);
			t1_old -> SetBranchAddress("beta_err",&LinkedLineStruct.beta_err);
			t1_old -> SetBranchAddress("ampl",&LinkedLineStruct.ampl);
			t1_old -> SetBranchAddress("ampl_err",&LinkedLineStruct.ampl_err);
			t1_old -> SetBranchAddress("tail",&LinkedLineStruct.tail);
			t1_old -> SetBranchAddress("tail_err",&LinkedLineStruct.tail_err);
			//t1_old -> SetBranchAddress("ratio",&LinkedLineStruct.ratio);
			//t1_old -> SetBranchAddress("ratio_err",&LinkedLineStruct.ratio_err);
			t1_old -> SetBranchAddress("step",&LinkedLineStruct.step);
			t1_old -> SetBranchAddress("step_err",&LinkedLineStruct.step_err);
			t1_old -> SetBranchAddress("cost",&LinkedLineStruct.cost);
			t1_old -> SetBranchAddress("cost_err",&LinkedLineStruct.cost_err);
		
			if( t1_old->GetBranch("chi2") ) t1_old -> SetBranchAddress("chi2", &LinkedLineStruct.chi2);
			if( t1_old->GetBranch("chi2ndof") ) t1_old -> SetBranchAddress("chi2ndof", &LinkedLineStruct.chi2ndof);
		
			if( t1_old->GetBranch("p_value") ) t1_old -> SetBranchAddress("p_value",&LinkedLineStruct.p_value);
			if( t1_old->GetBranch("p_value_ndof") ) t1_old -> SetBranchAddress("p_value_ndof",&LinkedLineStruct.p_value_ndof);
		
			if( t1_old->GetBranch("histo") ) t1_old -> SetBranchAddress("histo", &LinkedLineStruct.histo);
			if( t1_old->GetBranch("fit") ) t1_old -> SetBranchAddress("fit", &LinkedLineStruct.fit);
		
			for(int iLine=0; iLine<t1_old->GetEntries(); iLine++){
				t1_old->GetEntry(iLine);
				CalibLine tmpline = LinkedLineStruct;
				OldCalibLinesVector.push_back(tmpline);
			}
		}
	}
	
	
	vector<CalibLine> SavedCalibLineVector;
	
	for(int iLine=0; iLine<nlines; iLine++){
		string ans("");
		while(true){
			cout << "\nDo you want to fit the " << CalibLinesVec[iLine].massN << "-" << CalibLinesVec[iLine].element << " line at " << CalibLinesVec[iLine].litEn << " keV line? [y/n]\n" << "> ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("N") || ans!=string("n") || ans!=string("N")) break;
		}
		
		if(ans==string("y")  || ans==string("Y")){
			if(doFitBAT(MCAhisto,CalibLinesVec[iLine])){
				
				if(recreate){
					LinkedLineStruct = CalibLinesVec[iLine];
					t1->Fill();
					t1->AutoSave();
				}else{
					SavedCalibLineVector.push_back( CalibLinesVec[iLine] );
				}
				
			}
		}
	}
	
	
	if(!recreate){
		stringstream ss_tmp;
		//Merge the old and the new list of the calib line
		for(unsigned iLine=0; iLine<OldCalibLinesVector.size(); iLine++){
			ss_tmp.str(""); ss_tmp << OldCalibLinesVector.at(iLine).massN << OldCalibLinesVector.at(iLine).element << ((int)(OldCalibLinesVector.at(iLine).litEn+0.5));
			
			string oldname = ss_tmp.str();
			
			bool updated = false;
			for(unsigned kLine=0; kLine<SavedCalibLineVector.size(); kLine++){
				ss_tmp.str(""); ss_tmp << SavedCalibLineVector.at(kLine).massN << SavedCalibLineVector.at(kLine).element << ((int)(SavedCalibLineVector.at(kLine).litEn+0.5));
				
				string newname = ss_tmp.str();
				
				if(newname==oldname) updated = true;
			}
			
			if(!updated) SavedCalibLineVector.push_back( OldCalibLinesVector.at(iLine) );
			
		}
		
		for(unsigned iLine=0; iLine<SavedCalibLineVector.size(); iLine++){
			LinkedLineStruct = SavedCalibLineVector.at(iLine);
			t1->Fill();
		}
		t1->AutoSave();
	}
	
	if(outfile){
		outfile->Close();
		delete outfile;
	}
	
	return 0;
	
}
//===========================================================================//



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
	
	CalibLine LinkedLineStruct;//This is to address the trees.
	
	
	t1 -> Branch("massNum","string",&LinkedLineStruct.massN);
	t1 -> Branch("element","string",&LinkedLineStruct.element);
	t1 -> Branch("litEn",&LinkedLineStruct.litEn,"litEn/D");
	t1 -> Branch("litEn_err",&LinkedLineStruct.litEn_err,"litEn_err/D");
	t1 -> Branch("mean",&LinkedLineStruct.mean,"mean/D");
	t1 -> Branch("mean_err",&LinkedLineStruct.mean_err,"mean_err/D");
	t1 -> Branch("sigma",&LinkedLineStruct.sigma,"sigma/D");
	t1 -> Branch("sigma_err",&LinkedLineStruct.sigma_err,"sigma_err/D");
	t1 -> Branch("beta",&LinkedLineStruct.beta,"beta/D");
	t1 -> Branch("beta_err",&LinkedLineStruct.beta_err,"beta_err/D");
	t1 -> Branch("ampl",&LinkedLineStruct.ampl,"ampl/D");
	t1 -> Branch("ampl_err",&LinkedLineStruct.ampl_err,"ampl_err/D");
	t1 -> Branch("tail",&LinkedLineStruct.tail,"tail/D");
	t1 -> Branch("tail_err",&LinkedLineStruct.tail_err,"tail_err/D");
	//t1 -> Branch("ratio",&LinkedLineStruct.ratio,"ratio/D");
	//t1 -> Branch("ratio_err",&LinkedLineStruct.ratio_err,"ratio_err/D");
	t1 -> Branch("step",&LinkedLineStruct.step,"step/D");
	t1 -> Branch("step_err",&LinkedLineStruct.step_err,"step_err/D");
	t1 -> Branch("cost",&LinkedLineStruct.cost,"cost/D");
	t1 -> Branch("cost_err",&LinkedLineStruct.cost_err,"cost_err/D");
	
	t1 -> Branch("chi2",&LinkedLineStruct.chi2,"chi2/D");
	t1 -> Branch("chi2ndof",&LinkedLineStruct.chi2ndof,"chi2ndof/D");
	
	t1 -> Branch("p_value",&LinkedLineStruct.p_value,"p_value/D");
	t1 -> Branch("p_value_ndof",&LinkedLineStruct.p_value_ndof,"p_value_ndof/D");
	
	t1 -> Branch("histo", "TH1D", &LinkedLineStruct.histo);
	t1 -> Branch("fit", "TF1", &LinkedLineStruct.fit);
	
	if(!recreate){
		t1_old -> SetBranchAddress("massNum",&LinkedLineStruct.massN);
		t1_old -> SetBranchAddress("element",&LinkedLineStruct.element);
		t1_old -> SetBranchAddress("litEn",&LinkedLineStruct.litEn);
		t1_old -> SetBranchAddress("litEn_err",&LinkedLineStruct.litEn_err);
		t1_old -> SetBranchAddress("mean",&LinkedLineStruct.mean);
		t1_old -> SetBranchAddress("mean_err",&LinkedLineStruct.mean_err);
		t1_old -> SetBranchAddress("sigma",&LinkedLineStruct.sigma);
		t1_old -> SetBranchAddress("sigma_err",&LinkedLineStruct.sigma_err);
		t1_old -> SetBranchAddress("beta",&LinkedLineStruct.beta);
		t1_old -> SetBranchAddress("beta_err",&LinkedLineStruct.beta_err);
		t1_old -> SetBranchAddress("ampl",&LinkedLineStruct.ampl);
		t1_old -> SetBranchAddress("ampl_err",&LinkedLineStruct.ampl_err);
		t1_old -> SetBranchAddress("tail",&LinkedLineStruct.tail);
		t1_old -> SetBranchAddress("tail_err",&LinkedLineStruct.tail_err);
		//t1_old -> SetBranchAddress("ratio",&LinkedLineStruct.ratio);
		//t1_old -> SetBranchAddress("ratio_err",&LinkedLineStruct.ratio_err);
		t1_old -> SetBranchAddress("step",&LinkedLineStruct.step);
		t1_old -> SetBranchAddress("step_err",&LinkedLineStruct.step_err);
		t1_old -> SetBranchAddress("cost",&LinkedLineStruct.cost);
		t1_old -> SetBranchAddress("cost_err",&LinkedLineStruct.cost_err);
		
		if( t1_old->GetBranch("chi2") ) t1_old -> SetBranchAddress("chi2", &LinkedLineStruct.chi2);
		if( t1_old->GetBranch("chi2ndof") ) t1_old -> SetBranchAddress("chi2ndof", &LinkedLineStruct.chi2ndof);
		
		t1_old -> SetBranchAddress("p_value",&LinkedLineStruct.p_value);
		t1_old -> SetBranchAddress("p_value_ndof",&LinkedLineStruct.p_value_ndof);
		
		if( t1_old->GetBranch("histo") ) t1_old -> SetBranchAddress("histo", &LinkedLineStruct.histo);
		if( t1_old->GetBranch("fit") ) t1_old -> SetBranchAddress("fit", &LinkedLineStruct.fit);
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
				
				LinkedLineStruct = CalibLinesVec[iLine];
				
				t1->Fill();
				t1->AutoSave();
			}else{
				if(t1_old){
					for(int iEnt=0; iEnt<t1_old->GetEntries(); iEnt++){
						t1_old->GetEntry(iEnt);
						if( (LinkedLineStruct.massN==CalibLinesVec[iLine].massN) && (LinkedLineStruct.element==CalibLinesVec[iLine].element) && (LinkedLineStruct.litEn==CalibLinesVec[iLine].litEn)){
							t1->Fill();
							t1->AutoSave();
							break;
						}
					}
				}
			}
		}
	}
	
	
	if(outfile){
		if(t1) outfile->WriteTObject(t1);
		outfile->Close();
		delete outfile;
	}
	
	return 0;
	
}