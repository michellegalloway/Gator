#ifndef XURICH_LCELIKEL_HH
#define XURICH_LCELIKEL_HH

#include <vector>
#include <string>

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

#include "LikelihoodClass.hh"

#include "GatorStructs.h"

namespace Analysis{
	
	namespace GatorCalib{
	
		class GammaLineLikelihood: public Analysis::LikelihoodClass
		{

		public:
	
			GammaLineLikelihood();
			~GammaLineLikelihood(){;};
		
			
			//These methods resets the class to initial state. All the internal info are lost.
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
			void SetParStepSize(double step);
			void SetInitPars(vector<double> parVals);
			void SetMinuit2StepSize(double step){fMinuit2StepSize = step;};
		
			void MetropolisMLE();
		
			void Minuit2MLE(vector<double> parVal);
			void Minuit2MLE();
			
			//Miscellaneous
			//This is to scale back to ~0 the max log likelihood. Might be usefull for the numerical precision when used before than the minuit minimizer.
			void ScaleLogProbToMax(bool reset=false);
			
			
		private:
		
			int DefineParameters();
		
			vector<double> SelfParInit();
			
			
			//This computes the real likelihood must be implemented in a concrete class
			double LogProb(const vector<double>& par);
			
			//This return the opposite of the LogProb used in the Minuit2 minimizer
			static double MinLogProb(const double* parVal);
			
			
			void CalculatePValueLikelihood(vector<double> *pars);
			void CalculatePValueLikelihood(){CalculatePValueLikelihood(NULL);};
			
		
			TH1D *fHisto;
			int fNbins;
			
			CalibLine* fLine;
		
			//This is the scale of the step size for the MCMC proposals
			bool fExtParStep;
			double fParStepSize;
		
			//This is the data. Vector of doubles for each z slice.
			//fData[iSlice][iEntry]
			//The entries represent the ratio S1top/S1bot
			vector<vector<double>* >* fData;
			
			
			//This flag is true when the class is properly initialized for the MLE calculation
			bool fInit;
		
			ROOT::Minuit2::Minuit2Minimizer *fMinimizer;
			ROOT::Math::Functor *fMinFunctor;
		
			double fMinuit2StepSize;
		
			double fLogProbScaling;//This is always initialized to 0 and changed only by ScaleLobProbToMax() function 
		
			bool fMinimized;//It is true when at least one one of the MLE methods produced a good minimum
			double fPval, fPvalNDof; //Initialized to -1.0 in order to provide an additional flag that where not calculated
			
			bool fPvalues;//Flag to see if p-values have been calculated
			
		protected:
	
	
	
		};
		
	}//End of XurichLCE namespace
}//End of Analysis namespace





#endif
