#ifndef XURICH_LINELIKEL_HH
#define XURICH_LINELIKEL_HH

#include "LikelihoodClass.hh"
#include "GatorStructs.h"

#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TMath.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/DistFuncMathCore.h"
#include "Math/DistFuncMathMore.h"

#include <vector>
#include <string>



class GammaLineLikelihood: public LikelihoodClass
{
public:

	GammaLineLikelihood();
	~GammaLineLikelihood();

	
	//These methods reset the class to initial state. All the internal info are lost.
	int Init(TH1D* histo, CalibLine* line);
	
	//Getters
	double GetPval(vector<double> *pars);
	double GetPval(){ return GetPval(NULL); };
	
	double GetPvalNDof(vector<double> *pars);
	double GetPvalNDof(){ return GetPvalNDof(NULL); };
	
	int GetNDoF() const { return fNbins-fParsN; }
	
	double GetChi2(vector<double> *pars);
	double GetChi2() { return GetChi2(NULL); };
	
	double GetChi2NDof(vector<double> *pars);
	double GetChi2NDof() { return GetChi2NDof(NULL); };
	
	//Setters
	//void SetParStepSize(double step);
	void SetInitPars(vector<double> parVals);
	void SetMinuit2StepSize(double step){fMinuit2StepSize = step;};

	void MetropolisMLE();

	void Minuit2MLE(vector<double> &parVal);
	void Minuit2MLE();
	
	//Miscellaneous
	//This is to scale back to ~0 the max log likelihood. Might be usefull for the numerical precision when used before than the minuit minimizer.
	void ScaleLogProbToMax(bool reset=false);
	
	
private:
	
	virtual void MetropolisMLE(vector<double> *initpar);

	int DefineParameters();

	vector<double> SelfParInit();
	
	
	//This computes the real likelihood must be implemented in a concrete class
	double LogProb(const vector<double>& par);
	
	//This return the opposite of the LogProb used in the Minuit2 minimizer
	static double MinLogProb(const double* parVal);
	
	
	void CalculatePValueLikelihood(vector<double> *pars);
	void CalculatePValueLikelihood(){CalculatePValueLikelihood(NULL);};
	

	TH1D *fHisto;//Don't delete in the distructor
	int fNbins;
	
	CalibLine* fLine;//Don't delete in the distructor

	//This is the scale of the step size for the MCMC proposals
	bool fExtParStep;
	double fParStepSize;

	//This flag is true when the class is properly initialized for the MLE calculation
	bool fInit;

	ROOT::Minuit2::Minuit2Minimizer *fMinimizer;
	ROOT::Math::Functor *fMinFunctor;

	double fMinuit2StepSize;

	double fLogProbScaling;//This is always initialized to 0 and changed only by ScaleLobProbToMax() function 

	bool fMinimized;//It is true when at least one one of the MLE methods produced a good minimum
	double fPval, fPvalNDof; //Initialized to -1.0 in order to provide an additional flag that where not calculated
	
	bool fPvalues;//Flag to see if p-values have been calculated
	
	int fNsteps;
	
protected:
};


#endif
