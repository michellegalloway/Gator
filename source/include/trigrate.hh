#ifndef TRIGRATE_H
#define TRIGRATE_H


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TDatime.h>
#include <TMath.h>

//class SPEManager;

using namespace std;

class trigrate{

public:
	trigrate(); //Safe constructor
	trigrate(const Char_t* filename, const Char_t* workdir="./"); //With this contructor the file "triggerrate.txt" can be overwritten.
	virtual ~trigrate();
	
	TGraphErrors* MakeGraph(ULong_t startime, ULong_t endtime);
	TGraphErrors* MakeFullPlot();
	
	ULong_t GetTmin(){return m_minTstamp;};
	ULong_t GetTmax(){return m_maxTstamp;};
	Double_t GetEntry(Int_t i);
	Double_t GetEntryErr(Int_t i);
	Int_t Size(){return m_datasize;};
	
	
	//Sort all the 3 data vectors based on the timestamp
	void SortData();
	
	
	//Routine to write data in the static archive in append and in overwrite mode
	void WriteData(const Char_t* filename, Bool_t append);
	
	
	//Routine to print the last "num" dates
	void PrintDates(Int_t num);
	
	
	//This is the safe routine to be used to update the "triggerrate.txt" file. But can be also used in "non updating" mode. It saves the data of new found SPEs in the m_new... vectors.
	void FindSPE(const Char_t* dirname, Bool_t update_flg=false, Bool_t onlysource_flg=false); //The original one;
	void FindSPE(const Char_t* dirname, double StartCh, double EndCh, Bool_t update_flg=false, Bool_t onlysource_flg=false); //The new one version to set the start and the end channel
	
	//The input flag is to get from the manager the list of Spe files (sources) instead of SPE files (samples). The routine can be used only by the SPEManager.
	//void UpdateTrigRate(Bool_t onlysource=false);
	
	
	//This routine overwrite the "triggerrate.txt" (or equivalent "static trigger archive") file. Should be used only if you realize that data are corrupted after a certain time (it cancels the data after this time). After this routine should be run the FindSPE(..) routine to be sure that all the SPE files are in. Be carefull!!!
	void RepairTxt(const Char_t* SPEdir);
	
	
	//Usual routine to extract info from SPE data file. Originally designed to extract only the spectrum histogram.
	TH1F* LoadSPE(const Char_t* filename, Double_t& time, ULong_t& unixtime, std::string* descr_str=0);
	
	
	//Set the pointer to an SPEManager, which usually use this routine to give the pointer to itself
	//void SetSPEManager(SPEManager* pSPEmanager){m_pSPEmanager = pSPEmanager;};
	
	
	void SetWorkDir(const Char_t* dir);
	
	
private:
	//This routine loads the file's "triggerrate.txt" data. It should be called only by the constructor (for safety)
	void LoadTxt();
	
	std::string m_filename; //The name of the txt file like "triggerrate.txt"
	std::vector<ULong_t>* m_ptimestamp;
	std::vector<Double_t>* m_ptrigrate;
	std::vector<Double_t>* m_ptrigrate_err;
	
	//This are used for the update and for the repairing routine
	std::vector<ULong_t>* m_newtimestamp;
	std::vector<Double_t>* m_newtrigrate;
	std::vector<Double_t>* m_newtrigrate_err;
	
	const TDatime ch_date1; //Date when DAQ changed. -> From 12300-11 chs (11 are broken) to 16384 chs.
	
//	TGraphErrors* m_pgraph;
	
	string m_workdir; //The directory used as output dir
	
	Bool_t m_sorted;
	Bool_t m_txtfileflg; //This flag is true when the file's "triggerrate.txt" data are loaded. Can be changed only by LoadTxt() routine.
	
	Int_t m_datasize;
	
	ULong_t m_minTstamp, m_maxTstamp;
	
	//SPEManager* m_pSPEmanager;
};

#endif
