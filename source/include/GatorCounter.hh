#ifndef GATOR_COUNTER_HH
#define GATOR_COUNTER_HH

#include <string>

#include "GatorStructs.h"

class GatorCounter{

public:
	GatorCounter();
	~GatorCounter(){;};
	
	void SetBackgroundFlag(bool flag){mb_BGflag=flag;};
	void LineCounts(LineStruct &line);
	
	vector<LineStruct>* GetLines(){return mv_lines;};
	LineStruct* GetLine(int iLine);
	int GetNLines();
	
	
protected:
	
//In this two functions the two control regions are taken at the left and at the right of the central region, all with the same widths
	
	void ResetCounts();
	
	void SetPoissonCounting(){mb_nonIntCounts=false;};
	
	//Here I take into account when the edge of a region falls in the middle of a bin (not good to set up methods based on Poisson counting... Marc's way)
	void NonIntRegionCounts(LineStruct &line);
	
	//Here I take into account only integer counts (and the edges of the regions correspond always with bins boundaries)
	void PoissonRegionCounts(LineStruct &line);
	
	
	//If we are working or not whith the background
	bool mb_BGflag;
	
	//This flag is to consider non integer counts (when the region edge falls in the middle of a bin)
	bool mb_nonIntCounts; //[default: true]
	
	void SetDataSpect(TH1D* MCAspec){mp_sampleSpect=MCAspec;};
	
	void SetBackgroundSpect(TH1D* MCAspec){mp_backgroundSpect=MCAspec;};
	
	vector<LineStruct> *mv_lines;//Created by default as an empty vector
	
	TH1D *mp_sampleSpect; //Spectrum in MCA and not rebinned or rescaled
	TH1D *mp_backgroundSpect; //Spectrum in MCA and not rebinned or rescaled
	
};



#endif