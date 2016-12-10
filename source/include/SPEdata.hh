#ifndef SPEDATA_H
#define SPEDATA_H

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TH1D.h>
#include <TDatime.h>
#include <TMath.h>

//This class should just store the information contained in an SPE file. Not designed to do much more.
class SPEdata{

public:
	
	SPEdata(); //Can be used just to store data from other 
	SPEdata(string SPEfile); //This loads data from an SPE file
	
	~SPEdata();
	
	vector<double>* GetData(){return mv_pDataBins;};
	ULong_t GetStartUT(){return m_startunixtime;};
	double GetAcqTime(){return m_acqtime;};
	uint Size(){return m_nbins;};
	string GetDescription(){return ms_description;};
	
	TH1D* GetHisto(const string name=string(""),const string title=string(""));
	
	
	void SetData(vector<double>* databins); //This work only if data has not already been loaded from an SPE file. This to avoid the risk to screw up everything.
	void SetStartUT(ULong_t startunixtime){m_startunixtime=startunixtime;};
//	void SetSize(uint nbins){m_nbins=nbins;};
	void SetAcqTime(double acqtime){m_acqtime=acqtime;};
	void SetDescription(string descr){ms_descr=descr;};
	
private:

	//Routine to extract info from SPE data file. Used only in the constructor.
	void LoadSPEfile(string file);
	
	ULong_t m_startunixtime;
	double m_acqtime;
	vector<double>* mv_pDataBins; //Here is the data! Could create keep data in an TH1* but with this choice the class is much lighter (can create vector out of it.
	string ms_descr;
	uint m_nbins;
	
	bool flag_datafromfile; //Set on true only when the constructor that reads from the file is used.
	bool flag_dataloaded;
	
};

#endif