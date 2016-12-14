#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMinuit.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "GatorCalibFitters.h"

#include "GammaLineLikelihood.hh"

using namespace std;
using Analysis::Param;

namespace Analysis{
	namespace GatorCalib{
		
		
		GammaLineLikelihood::GammaLineLikelihood():
		{
			fParsN = 7;
			fParams = vector<Param*>(fParsN);
			fParSteps = new vector<double>(fParsN);
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				fParams.at(iPar)=NULL;
			}
			fCalibLine = NULL;
			fInit =false;
			fExtParStep=false;
			fParStepSize=0.01;
			fHisto=NULL;
			fLine=NULL;
			fMinimizer=NULL;
			fMinFunctor=NULL;
			fMinuit2StepSize=0.001;
			fLogProbScaling=0.;
			fPvalues = false;
			fPval=-1.0;
			fPvalNDof=-1.0;
		}

		int GammaLineLikelihood::Init(TH1D* hist, const CalibLine* line)
		{
			//The class must be reset
			if(fInit){
				fInit=false;
			}
			fHisto=NULL;
			fLine=NULL;
			fExtParStep=false;
			fParStepSize=0.01;
			fHisto=NULL;
			fMinimizer=NULL;
			fMinFunctor=NULL;
			fMinuit2StepSize=0.001;
			fLogProbScaling=0.;
			fWriteOutput=false;
			fPval=-1.0;
			fPvalNDof=-1.0;
			
			
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				if( fParams.at(iPar) ){
					delete fParams.at(iPar);
					fParams.at(iPar)=NULL;
				}
			}
			
			int errcode = DefineParameters();
			if(errcode!=0){
				cerr << "GammaLineLikelihood::Init(...) --> ERROR: DefineParameters() failure!\nExiting with error code " << errcode << endl << endl;
				return(errcode);
			}else{
				fLine = line;
				fHisto = histo;
				fNbins = histo->GetNbinsX();
			}
			
			fInit=true;
			return(0);
	
		}


		int GammaLineLikelihood::DefineParameters()
		{
			cout << "\nEntering in GammaLineLikelihood::DefineParameters()" << endl;
			
			double Val, transfVal;
			
			//Here I define the parameters
			Param *tmp_par = new Analysis::Param("Mean");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->mean;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(0) = tmp_par;
			
			
			Param *tmp_par = new Analysis::Param("Ampl");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->ampl;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(1) = tmp_par;
			
			Param *tmp_par = new Analysis::Param("Tail");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->tail;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(2) = tmp_par;
			
			Param *tmp_par = new Analysis::Param("Sigma");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->sigma;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(3) = tmp_par;
			
			Param *tmp_par = new Analysis::Param("Beta");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->beta;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(4) = tmp_par;
			
			Param *tmp_par = new Analysis::Param("Step");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->step;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(5) = tmp_par;
			
			Param *tmp_par = new Analysis::Param("Const");
			tmp_par->SetLowerLimit(0);
			tmp_par->SetUpperLimit(1);
			//Transform the parameter from the (0,inf) range to the (0,1) range
			val = line->cost;
			transfVal = val/(1+val);
			tmp_par->SetValue(transfVal);
			fParams.at(6) = tmp_par;
			
			if(fTrackVals) tmp_par->TrackVals();
	
			fMaxPar = vector<double>(fParsN);
			if(fVerbosity>=1){
				cout << "-> " << fParsN << " parameters defined and initialised." << endl;
			}
			
			fParInit = true;
			
			cout << "Exiting from GammaLineLikelihood::DefineParameters()\n" << endl;
			
			return(0);
		}


		void GammaLineLikelihood::SetParStepSize(double step)
		{
			fParStepSize=step;
			if(!fParSteps){
				fParSteps = new vector<double>(fParsN);
			}
			
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				fParSteps->at(iPar) = fParStepSize;
			}
			
			fExtParStep=true;
		}


		void GammaLineLikelihood::SetInitPars(vector<double> parVal)
		{
			vector<double> *initpar = new vector<double>( parVal.begin(), parVal.begin()+parVal.size() );
			Analysis::LikelihoodClass::SetInitPars(initpar);
			return;
		}


		vector<double> GammaLineLikelihood::SelfParInit()
		{
			//This function is not needed as the parameters are already initialized where they are defined
			vector<double> par(fParsN);
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				par.at(iPar) = fParams.at(iPar)->GetValue();
			}
			return par;
		}


		double GammaLineLikelihood::LogProb(const vector<double>& par)
		{
			//Vector to transform back the parameters from (0,1) interval to (0,infinity)
			static double *truepars = new double[fParsN];
			
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				truepars[iPar] = par.at(iPar)/(1.-par.at(iPar));
			}
			
			//Calculate the log-likelihood
			if(fEngineVerb>=4){
				cout << "\n====> Entering in GammaLineLikelihood::LogProb(...)." << endl;
			}
			double logprob =0;
	
			for(int iBin=1; iBin<=fNbins; iBin++){
				const double xmin = fHisto->GetBinLowEdge(iBin);
				const double xmax = fHisto->GetBinLowEdge(iBin+1);
				
				//Interpolate the function in the central interval and compare with the bin contents
				const double expCounts = (peakFitFunc(&xmax, truepars) + peakFitFunc(&xmin, truepars))*(xmax-xmin)/2.;
				
				const double obsCounts = fHisto->GetBinContent(iBin);
				
				logprob += obsCounts*log(expCounts) - expCounts;
				
			}
			
			
			if(fEngineVerb>=4){
				cout << "====> Exiting from GammaLineLikelihood::LogProb(...) with value: " << logprob << endl;
			}
			return logprob-fLogProbScaling;
		}
		
		
		void GammaLineLikelihood::ScaleLogProbToMax(bool reset)
		{
			if(reset){
				fLogProbScaling = 0;
				return;
			}
	
			fLogProbScaling = fMaxLogProb;
			return;
		}


		double GammaLineLikelihood::MinLogProb(const vector<double> *parVal)
		{
			//This return the opposite of the LogProb used in the Minuit2 minimizer
	
			Analysis::GatorCalib::GammaLineLikelihood* this_ptr = static_cast<Analysis::GatorCalib::GammaLineLikelihood*>(Analysis::LikelihoodClassHolder::instance());
	
	
			//vector<double> par = *parVal;
			return -(this_ptr->LogProb(*parVal));
		}


		void GammaLineLikelihood::MetropolisMLE()
		{
			if(!fInit) return;
	
			Analysis::LikelihoodClass::MetropolisMLE();
			fMinimized = true;
			return;
		}


		void GammaLineLikelihood::Minuit2MLE(vector<double> parVal)
		{
			cout << "Analysis::GammaLineLikelihood::Minuit2MLE() function" << endl;
			
			if(fMinimizer){
				delete fMinimizer;
				fMinimizer=NULL;
			}
			if(fMinFunctor){
				delete fMinFunctor;
				fMinFunctor = NULL;
			}
			
			fMinimizer = new ROOT::Minuit2::Minuit2Minimizer();
			
			// set tolerance , etc...
			fMinimizer->SetMaxFunctionCalls(100000); // for Minuit/Minuit2 
			fMinimizer->SetTolerance(0.001);
			
			//Setverbose in testing phase
			if(fEngineVerb<=1){
				fMinimizer -> SetPrintLevel(-1);
			}
			else if(fEngineVerb==2){
				fMinimizer -> SetPrintLevel(0);
			}
			else if(fEngineVerb>2){
				fMinimizer -> SetPrintLevel(1);
			}
			
			
			//Set the instance before using the the static method MinLogProb in the minimizer
			Analysis::LikelihoodClassHolder::instance(this);
			//fMinFunctor = new ROOT::Math::Functor(this,&XurichAnalysis::XuLCElikelihood::MinLogProb,1);
			fMinFunctor = new ROOT::Math::Functor(&Analysis::GatorCalib::GammaLineLikelihood::MinLogProb, fParsN);
			
			fMinimizer->SetFunction(*fMinFunctor);
			
			
			vector<double> truepars(fParsN);
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				//Transform back the parameters from (0,1) interval to (0,infinity)
				double trueparval = par.at(iPar)/(1.-par.at(iPar));
				fMinimizer->SetVariable(iPar, fParams.at(iPar)->GetName(), trueparval, fMinuit2StepSize);
			}
			
			
			if(fEngineVerb>=2) cout << "=>Start the Minuit2 minimization." << endl;
			
			if(!fMinimizer->Minimize()){
				cout << "WARNING: Minimization in Minuit 2 failed with error code " << fMinimizer->Status() << endl;
				return;
			}
			
			if(fEngineVerb>=2) cout << "=>Start the Minuit2 hessian calculation." << endl;
			
			if(!fMinimizer->Hesse()){
				cout << "WARNING: Hessian matrix in Minuit 2 failed with error code " << fMinimizer->Status() << endl;
				return;
			}
			
			cout << "\nMinuit 2 found the minimum of parameter <" << fParams.at(0)->GetName() << ">:"  << endl;
			cout << "Min value = " << fMinimizer->MinValue() << endl;
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				cout << fParams.at(iPar)->GetName() << " MLE = (" << (fMinimizer->X())[iPar] << " +- " << sqrt(fMinimizer->CovMatrix(iPar,iPar)) << ")" << endl;
			}
			
			
			line->mean = (fMinimizer->X())[0];
			line->mean_err = sqrt(fMinimizer->CovMatrix(0,0));
			
			line->ampl = (fMinimizer->X())[1];
			line->ampl_err = sqrt(fMinimizer->CovMatrix(1,1));
			
			line->tail = (fMinimizer->X())[2];
			line->tail_err = sqrt(fMinimizer->CovMatrix(2,2));
			
			line->sigma = (fMinimizer->X())[3];
			line->sigma_err = sqrt(fMinimizer->CovMatrix(3,3));
			
			line->beta = (fMinimizer->X())[4];
			line->beta_err = sqrt(fMinimizer->CovMatrix(4,4));
			
			line->step = (fMinimizer->X())[5];
			line->step_err = sqrt(fMinimizer->CovMatrix(5,5));
			
			line->cost = (fMinimizer->X())[6];
			line->cost_err = sqrt(fMinimizer->CovMatrix(6,6));
			
			line->chi2 = GetChi2();
			line->chi2ndof = GetChi2NDof();
			
			line->p_value = GetPValue();
			line->p_value_ndof = GetPValueNDoF();
			//cout << "Variance of the parameter = " << fMinimizer->CovMatrix(0,0) << ". Root square = " << sqrt(fMinimizer->CovMatrix(0,0)) << endl;
			//cout << "Error = " << (fMinimizer->Errors())[0] << endl;
			
			/*
			if(fWriteOutput){
				fResultsOut.open("LCEResults.txt");
				fResultsOut << "\nMinuit 2 found the minimum of parameter <" << fParams.at(0)->GetName() << ">:"  << endl;
				fResultsOut << "MLE = " << (fMinimizer->X())[0] << endl;
				fResultsOut << "Min value = " << fMinimizer->MinValue() << endl;
				fResultsOut << "Variance of the parameter = " << fMinimizer->CovMatrix(0,0) << ". Root square = " << sqrt(fMinimizer->CovMatrix(0,0)) << endl;
				fResultsOut << "Error = " << (fMinimizer->Errors())[0] << endl;
				fResultsOut.close();
			}*/
			
			cout << "Exiting from Analysis::GammaLineLikelihood::Minuit2MLE(...) function." << endl;
			
			return;
		}

		
		void GammaLineLikelihood::CalculatePValueLikelihood(vector<double> *pars)
		{
			
			if(!fMinimized) return;
			
			vector<double> truepars(fParsN);
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				truepars.at(iPar) = fMaxPar.at(iPar)/(1.-fMaxPar.at(iPar));
			}
			
			// initialize test statistic -2*lambda
			double logLambda = 0.0;

			for(int iBin = 1; iBin <= fNbins; ++iBin) {
				// get the number of observed events
				const double y = fHisto->GetBinContent(ibin);
				
				// get the number of expected events using the interpolation of the function in the bin "iBin"
				const double xMin = fHisto->GetBinLowEdge(iBin);
				const double xMax = fHisto->GetBinLowEdge(iBin+1);
				
				const double yexp = ( peakFitFunc(&xMax, &(truepars.at(0))) + peakFitFunc(&xMin, &(truepars.at(0))) )*(xMax-xMin)/2.;
				
				// get the contribution from this datapoint
				if (y == 0)
					logLambda += yexp;
				else
					logLambda += yexp - y + y * log(y / yexp);
			}

			// rescale
			logLambda *= 2.0;

			//p value from chi^2 distribution, returns zero if logLambda < 0
			fPValue = TMath::Prob( logLambda, fNbins );
			fPValueNDoF = TMath::Prob( logLambda, GetNDoF() );
			
			fPvalues = true;
			
			// no error
			return;
		}
		
		
		double GammaLineLikelihood::GetPval(vector<double> *pars)
		{
			if(pars){
				if(!pars->size()==fParsN) return -1.0;
				CalculatePValueLikelihood(pars);
				return fPValue;
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if(!fPvalues) return -1.0;
			
			return fPValue;
		}
		
		
		double GammaLineLikelihood::GetPvalNDof(vector<double> *pars)
		{
			if(pars){
				if(!pars->size()==fParsN) return -1.0;
				CalculatePValueLikelihood(pars);
				return fPValueNDoF;
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if(!fPvalues) return -1.0;
			
			return fPValueNDoF;
		}
		
		
		double GetChi2(vector<double> *pars)
		{
			if(pars){
				if(!pars->size()==fParsN) return -1.0;
				CalculatePValueLikelihood(pars);
				if( !( (fPValueNDoF>0.0) && (<1.0) ) ) return -1.0;
				return ROOT::Math::chisquared_quantile_c( fPValueNDoF, (double)GetNDoF() );
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if( !( (fPValueNDoF>0.0) && (<1.0) ) ) return -1.0;
			
			return ROOT::Math::chisquared_quantile_c( fPValueNDoF, (double)GetNDoF() );
		}
		
		
		double GetChi2NDof(vector<double> *pars) const
		{ 
			if(pars){
				if(!pars->size()==fParsN) return -1.0;
				double chi2 = GetChi2(pars);
				if(chi2<=0.) return -1.0;
				return chi2/GetNDoF();
			}
			
			double chi2 = GetChi2();
			if(chi2<=0.) return -1.0;
			return chi2/GetNDoF();
		}
	}
}



