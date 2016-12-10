#ifndef GATOR_PARAM_HH
#define GATOR_PARAM_HH

#include <string>
#include <vector>
#include "TH1D.h"

using namespace std;


class GatorParam{

public:
	GatorParam(string name, double upEdge, bool trackvals=false);
	GatorParam(string name, double lowEdge, double upEdge, bool trackvals=false);
	~GatorParam(){;};
	
	//Getters
	const string & GetName() const { return fName; };
	inline double GetUpperLimit() const { return fUpperLimit; };
	inline double GetLowerLimit() const { return fLowerLimit; };
	//bool FillHistogram(); //This must be implemented later ()
	inline double GetValue() const { return fValue; };
	inline double GetLogProb() const { return fValue; };
	unsigned GetNbins() const { return fNbins; };
	TH1D* GetHisto() const { return fHist;};
	vector<vector<double> > GetTrackingVals() const { return fValsChain; };
	
	//Setters
	void SetName(string name) { fName = name; };
	bool SetValue(double value, double likelihood=0);
	void SetUpperLimit(double limit) { fUpperLimit = limit; };
	void SetLowerLimit(double limit) { fLowerLimit = limit; };
	void SetLimits(double lowerlimit, double upperlimit);
	void SetHist(bool flag=true) { fHistoFlag=flag; };
	void SetNbins(unsigned nbins) { fNbins = nbins; };
	void ForceValue(double value, double logprob);
	
	//Misc
	inline bool FillHisto() const { return fHistoFlag; };
	inline bool IsFixed() const { return fFixed; };
	inline bool IsTracking() const { return fTrackVals; };
	void TrackVals(bool trackval=true) { fTrackVals=trackval; };
	void Fix(double value);
	void Unfix(){ fFixed = false; };
	//void grams(bool flag) { fFillHistograms = flag; };
	bool IsAtLimit(double value) const { return ((fLowerLimit == value) || (value == fUpperLimit)); } ;
	int IsOutOfRange(double value) const;
	
private:
	
	string fName;
	
	double fUpperLimit;
	double fLowerLimit;//This should always be zero
	bool fFixed; //If the value is fixed
	double fValue; //The current value of the parameter
	double fLogProb;
	bool fHistoFlag; //If the histogram is filled
	unsigned fNbins; //Number of bins for the histogram
	bool fTrackVals;
	bool fLockLowEdge;//If you want the lower edge locked to zero (for bayesian for example)
	
	//Here I could dump the profiled likelihood or the marginal posterior pdf for the parameter.
	TH1D* fHist;
	//In the first vector are saved the values of the parameter and in the second the values of the likelihood (or the log-likelihood or the posterior pdf values)
	vector<vector<double> > fValsChain;
	
	
	
};


#endif