#ifndef GATOR_TREEDATALOADER_HH
#define GATOR_TREEDATALOADER_HH

#include <vector>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"

#include "GatorParam.hh"
#include "GatorCounter.hh"

//This class is made to load from a tree the spectra saved in arrays. Devoloped mostly for analysis of simulated data.
class GatorTreeDataLoader: public virtual GatorCounter{

public:
	//This is the real constructor
	GatorTreeDataLoader( string sourcefilename, string treename=string("t1") );
	
	//This constructor is to keep this class inactive
	GatorTreeDataLoader();
	
	inline bool IsTreeDataLoaded() const { return fDataLoaded; };
	
	
	//Getters
	inline vector<double>* GetLinesBgLevels(){ return &fBgLevel; };
	inline double GetSignalLev(){ return fSigLevel; };
	inline vector<double>* GetLinesSigma(){ return &fSigma; };
	inline int GetTreeEntries(){ return fTreeEntries; };
	
	~GatorTreeDataLoader(){;};
	
	
private:
	
	//This function makes the parameters out of line numbers. Should be used in the constructor. Returns an error code
	int LoadTree();
	bool fDataLoaded;
	//This computes the real likelihood
	//double LogProb(vector<double>& par);
	
	//This is to be used in the EM algorithm
	//double LogProbComplete(vector<double>& par); 
	
	//This is the method that returns the new proposed parameters
	//vector<double> MetrProp(vector<double>& par);


protected:
	unsigned fSize; //Is the lenght of simulated lines spectrum (read from the root file)
	
	string fSourceFileName;
	string fTreeName;
	//Here are stored ALL the parameters
	//vector<GatorParam*> fParams;
	unsigned m_nLines;
	//double fWidthScale; //This is to define the windows width respect to the line's width [default = 3.0]
	//double fAbsPrec; //Userfull Only when I use the EM algorithm
	//unsigned fMaxIter; //Max number of iterations to perform before returnig a result from Likelihood maximization
	//Nominal levels
	double fSigLevel; //Expressed in counts (number of decays)
	double fSnrLevel; //Expressed as the (SigLev/Sigma)/BgLevel (from the first of the gamma lines)
	vector<double> fSigma;//Expressed in bins
	vector<double> fBgLevel; //Expressed in counts per bin
	vector<double> fEffic;
	//bool fMultipleLines;
	
	//unsigned m_nPar;
	int fTreeEntries; //The number of the tree's entries
	bool fTreeDataLoaderActive; //If this class is active
	bool fSingleLine;
};

#endif