#ifndef GATOR_SAMPLEANALYSIS_HH
#define GATOR_SAMPLEANALYSIS_HH

#include <string>

#include "TCanvas.h"
#include "TPad.h"
#include "TH2F.h"


#include "GatorStructs.h"
#include "GatorCounter.hh"
#include "GatorDataLoader.hh"



TH2F* makeGraphSupport(double fMinimum1, double fMaximum1, double fMinimum2, double fMaximum2);


class GatorSampleAnalysis: public GatorDataLoader{

public:
	
	GatorSampleAnalysis(string workdir);
	~GatorSampleAnalysis(){;};
	
	TH1F* GetSampleEnSpect(){return mp_sampleEnSpect;};
	TH1F* GetBgEnSpect(){return mp_backgroundEnSpect;};
	
	TGraphErrors* GetDetectPlot(){return gr_DetActiv;};
	TGraphErrors* GetUpLimPlot(){return gr_DetActiv;};
	
	
	//This is to use the background even if it is not detected
	void ForcePeakBg(bool pkbg=true){mb_BgForced=pkbg;};
	
	void SetQuantityUnit(string unit){m_unit=unit;};
	
	//This can be the amount of sample (nPMTs for example)
	void SetSampleQuntity(double quant){md_quant=quant;};
	
	void SetLogFile(string logffile){m_logfilename=logffile;};
	
	void DrawActivityPlots(TPad* pad=NULL);
	
	void WriteOutputTable(string txtfile);
	
	void WriteEffTable(string txtfile);
	
	//Put a tree and the spectra inside the root file
	void WriteRootFile(string filename, string mode=string("update"));
	
private:
	
	//This fills the proper fields in the LineStruct (after some controls)
	void doAnalysis();
	
	string m_workdir; //The output directory
	string m_logfilename;
	
	string m_unit;
	
	double md_quant;
	
	TH1F *mp_sampleEnSpect; //Spectrum in MCA and not rebinned or rescaled
	TH1F *mp_backgroundEnSpect; //Spectrum in MCA and not rebinned or rescaled
	
	TGraphErrors *gr_DetActiv;
	TGraphErrors *gr_ULActiv;
	
	bool mb_dataAnalysed;
	
	bool mb_BgForced;//This is to use the background even if it is not detected
};

#endif