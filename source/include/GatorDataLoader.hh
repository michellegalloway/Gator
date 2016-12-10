#ifndef GATOR_DATALOADER_HH
#define GATOR_DATALOADER_HH

#include <string>

#include "GatorStructs.h"
#include "GatorCounter.hh"
#include "screenfncs.h"


class GatorDataLoader: public virtual GatorCounter{

public:
	
	GatorDataLoader();
	~GatorDataLoader(){;};
	
	void SetSampleDataDir(string dir);
	void SetBackgroundDataDir(string dir);
	void SetWidthScale(double scale){md_WidthScale=scale;};
	
	virtual void SetLeftCompRegion(double lowEn, double upEnr, int iLine);
	virtual void SetSignalRegion(double lowEn, double upEnr, int iLine);
	virtual void SetRightCompRegion(double lowEn, double upEnr, int iLine);
	
	virtual bool LoadLines(string listname);
	
	virtual bool LoadData(string dir=string(""));
	bool LoadBackground(string dir=string(""));
	
	bool SetSampleCalib(string dir);
	bool SetBackgroundCalib(string dir);
	
	virtual void DefaultIntervals();
	

protected:
	
	//map<string, int> linesmap; //This is to call lines by name instead from their number
	
	virtual void DefaultIntervals(LineStruct &line);
	
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