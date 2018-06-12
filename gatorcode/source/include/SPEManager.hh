#ifndef SPEMANAGER
#define SPEMANAGER

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TMath.h"

#include "trigrate.hh"


using namespace std;

class SPEManager{


public:

	SPEManager(string trigratename="triggerrate.txt", string workdir="", string spelistname="spedynlist");
	virtual ~SPEManager();
	
	//This is the script that should be run from the Gator slow control. It can be run in loop mode but before of that it must be debugged carefully.
	void runSPEManager();
	
	
	TH1F* GetADCSpectrum();
	
	
	//This gives the spectrum only rescaled by daily time -> this can be rescaled by bin width after rebinning in order to get on y-axis cnts/(sec*keV)
	TH1F* GetEnergySpectrum();
	
	
	//General getters
	int GetNSpectra(){return num_of_spectra;};
	Double_t GetAqTime(){return aqtime;};
	ULong_t GetLastTimestamp(){return last_timestamp;};
	string GetStaticArchive() {return cm_trigratename;};
	string GetWorkDir() {return cm_workdir;};
	vector<string>* GetNewSPEVector() {return &newspefiles;};
	vector<string>* GetNewSourceVector() {return &newsourcefiles;};
	string GetSampleDescr(){return m_sampledescr;};
	
	//Mandatory Setters
	void SetSPEdir(string dir);
	void SetSourcedir(string dir);
	

	
private:
	
	Bool_t FindNewSample();
	void FindNewSource(); //Not usefull because the Spe files (sources) don't go in the "spe_reports" directory
	void TrigrateUpdate();
	
	Bool_t CheckDescription(vector<string>& newspefiles);
	
	void ADCspectrumBuild();
	
	void EnergySpectrumBuild();
	
	
	string cm_trigratename; //The static list "triggerrate.txt"
	string cm_spelistname; //Name of the dynamic list with the SPE file names that build up the plot
	string cm_workdir; //The directory of work
	string cm_SPEdir; //The name of the dir where all the SPE files are moved
	string cm_Sourcedir; //The name of the dir where all the Spe (sources) files are moved
	const TDatime cm_date1; //Date when DAQ changed. -> From 12300-11 chs (11 are broken) to 16384 chs.
	trigrate* pm_trigrateobj;
	
	vector<string> newspefiles; //This is the vector with the SPE filenames containing the last description
	vector<string> newsourcefiles;
	string m_newsampledescr;
	string m_sampledescr;
	
	TH1F* m_pADCspectrum;
	
	TH1F* m_pEnergySpectrum;
	
	int num_of_spectra;
	Double_t aqtime;
	ULong_t last_timestamp;
	
	Bool_t newsamp_flg;
	Bool_t SPEdir_flg;
	Bool_t Sourcedir_flg;
	
	
};

#endif
