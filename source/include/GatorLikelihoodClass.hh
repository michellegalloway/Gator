#ifndef GATOR_LIKELIHOODCLASS_HH
#define GATOR_LIKELIHOODCLASS_HH

#include <vector>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"

#include "GatorParam.hh"
#include "GatorCounter.hh"
#include "GatorDataLoader.hh"
#include "GatorTreeDataLoader.hh"

//This a container to store the results of the log-likelihood and of the parameter values
class GatorLikelihoodResults{

public:
	double logprob;
	vector<double> par;
};


//This class is made to load from a tree the spectra saved in arrays. Devoloped mostly for analysis of simulated data.
class GatorLikelihoodClass: public GatorTreeDataLoader, public GatorDataLoader{

public:
	
	enum GatorDataType{
		kTreeData,
		kSpeData
	};
	
	//This constructor is to be used when the spectra are saved in a TTree. Should be used in algorithms testing and in limits and sensitivity derivations
	GatorLikelihoodClass( string sourcefilename, string treename );
	//This constructor is used to load real data
	GatorLikelihoodClass( string SPEdir ){;};
	~GatorLikelihoodClass(){;};
	
	
	virtual void AnalyzeTree(){;};
	
	//Like the method before but only one entry of the tree. The file here is opened in red mode only and the trace won't be saved (usefull for debugging and testing!)
	virtual void AnalyzeTreeEntry(int entry){;};
	
	//From here I can take all the values for the LogProb given a fixed value for the signal (or snr)
	virtual vector<double> TreePLs(double signal){return vector<double>(0);};
	
	//Getters
	inline vector<GatorParam*>& GetParameters() { return fParams; };
	GatorParam* GetParameter(unsigned iPar);
	unsigned GetNParameters() const { return m_nPar;};
	inline double GetWidthScale(){ return fWidthScale; };
	virtual double GetMaxLProb(){ return fMaxLogProb; };
	
	virtual void SingleLineMLE(){;};
	virtual void MetropolisMLE(vector<double> *initpar=NULL){;};
	
	//Setters
	//inline void SetWidthScale(double scale){ fWidthScale=scale; };
	//void SetAbsPrec(double prec){ fAbsPrec=prec;};
	//void SetMaxSteps(unsigned maxSteps){ fMaxIter=maxSteps; };
	//void SetTrackVals(bool trackvals=true){ fTrackVals=trackvals; };
	//This is to set the width of the regions respect to the sigma of the line
	inline void SetWidthMult(double multip){ fWidthScale=multip; };
	
private:
	
	
	
	
	
protected:
	
	//This function initialize the regions of the lines. To be defined by the user.
	virtual int InitLines()=0;
	
	//This method is responsible to prepare the data, create a new tree and make a full analysis calling other methods.
	virtual int DefineParameters()=0;
	
	//This computes the real likelihood must be implemented in a concrete class
	virtual double LogProb(vector<double>& par)=0;
	
	GatorDataType fDataType;
	
	double fWidthScale; //This is to define the windows width respect to the line's width [default = 3.0]
	
	unsigned m_nPar;
	TRandom3 *fRndGen;
	
	//Maximization stuff
	TMinuit *fMinuit;
	vector<double> fMaxPar;
	double fMaxLogProb;
	
	double fMinuitArglist[2];
	
	bool fSingleLine; //When I have a single line I can exploit the analytical solutions for the MLE
	bool fLinesInit; //Flag to chech if the lines regions are properly initialized
	bool fParInit; //Flag to check if the parameters are properly initialized
	
	vector<GatorParam*> fParams; //Here are stored ALL the parameters
	//unsigned m_nLines;
	
	
	//This stuff is redundant, but it is confy to use in the calculators (for speed)
	//vector<double> *w;
	vector<double> *Lcnts;
	vector<double> *Gcnts;
	vector<double> *Rcnts;
	vector<double> *ParSteps;
	
};


//I need this in order to use properly Minuit
namespace {
	class GatorLikelihoodClassHolder{
	private:
		GatorLikelihoodClass * global_this;

	public:
		
		GatorLikelihoodClassHolder():
		    global_this(NULL)
			{
			}
			
		static GatorLikelihoodClass * instance(GatorLikelihoodClass * obj = NULL){
			static GatorLikelihoodClassHolder result;
			if (obj) result.global_this = obj;
			return result.global_this;
		}
	};
}


#endif