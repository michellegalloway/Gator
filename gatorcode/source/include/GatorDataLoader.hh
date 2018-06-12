#ifndef GATOR_DATALOADER_HH
#define GATOR_DATALOADER_HH

#include <string>

#include "GatorStructs.h"
#include "GatorCounter.hh"

TH1F* loadSpe(const char* dir, Double_t& aqtime);


class GatorDataLoader: public GatorCounter{

public:
	
	GatorDataLoader();
	~GatorDataLoader(){;};
	
	void SetSampleDataDir(string dir);
	void SetBackgroundDataDir(string dir);
	void SetWidthScale(double scale){md_WidthScale=scale;};
	
	void SetLeftCompRegion(double lowEn, double upEnr, int iLine);
	void SetSignalRegion(double lowEn, double upEnr, int iLine);
	void SetRightCompRegion(double lowEn, double upEnr, int iLine);
	
	bool LoadLines(string listname);
	
	bool LoadData(string dir=string(""));
	bool LoadBackground(string dir=string(""));
	
	bool SetSampleCalib(string dir);
	bool SetBackgroundCalib(string dir);
	
	void DefaultIntervals();
	
protected:
	
	//map<string, int> linesmap; //This is to call lines by name instead from their number
	
	void DefaultIntervals(LineStruct &line);
	
	double md_aqtime, md_Bgtime;
	
	double md_WidthScale;
	
	string m_datadir;
	string m_datacalibdir;
	TF1 *mp_data_calib_fw; //This to make the energy spectra
	TF1 *mp_data_calib_bw; //This to setup the lineMCA channel and its sigma
	TF1 *mp_data_MCA_resol; //Resolution for this dataset
	
	string m_Bgdir;
	string m_Bgcalibdir;
	TF1 *mp_Bgcalib_fw; //This to make the energy spectra
	TF1 *mp_Bgcalib_bw; //This to setup the lineMCA channel and its sigma
	TF1 *mp_BgMCA_resol; //Resolution for this dataset
	
	bool mb_dataloaded, mb_Bgloaded, mb_linesloaded;
};


#endif