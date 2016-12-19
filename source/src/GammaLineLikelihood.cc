#include <limits>
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
		
		
		GammaLineLikelihood::GammaLineLikelihood()
		{
			fParsN = 0;
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
			fNsteps = 100000;
			fMaxIter = 10000000;
		}
		
		
		GammaLineLikelihood::~GammaLineLikelihood()
		{
			if(fMinimizer) delete fMinimizer;
			if(fMinFunctor) delete fMinFunctor;
			
			//Protected members of the base class
			if(fRndGen) delete fRndGen;
			if(fMinuit) delete fMinuit;
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				if(fParams.at(iPar)) delete fParams.at(iPar);
			}
			if(fInitPars) delete fInitPars;
			if(fParSteps) delete fParSteps;
		}
		
		
		int GammaLineLikelihood::Init(TH1D* histo, CalibLine* line)
		{
			//The class must be reset
			if(fInit){
				fInit=false;
			}
			fHisto=NULL;
			fExtParStep=false;
			fParStepSize=0.01;
			fHisto=NULL;
			if(fMinimizer){
				delete fMinimizer;
				fMinimizer=NULL;
			}
			if(fMinFunctor){
				delete fMinFunctor;
				fMinFunctor=NULL;
			}
			fMinuit2StepSize=0.001;
			fLogProbScaling=0.;
			fPval=-1.0;
			fPvalNDof=-1.0;
			
			
			fLine = line;
			
			for(unsigned iPar=0; iPar<fParams.size(); iPar++){
				if( fParams.at(iPar) ){
					delete fParams.at(iPar);
				}
			}
			fParams.clear();
			
			int errcode = DefineParameters();
			if(errcode!=0){
				cerr << "GammaLineLikelihood::Init(...) --> ERROR: DefineParameters() failure!\nExiting with error code " << errcode << endl << endl;
				fLine = NULL;
				return(errcode);
			}else{
				fHisto = histo;
				fNbins = histo->GetNbinsX();
			}
			
			fInit=true;
			return(0);
	
		}


		int GammaLineLikelihood::DefineParameters()
		{
			cout << "\nEntering in GammaLineLikelihood::DefineParameters()" << endl;
			
			double val, transfVal;
			
			Param *tmp_par;
			
			//Here I define the parameters
			tmp_par = new Analysis::Param("Mean");
			tmp_par->SetLowerLimit(0);
			val = fLine->mean;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/1000);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Ampl");
			tmp_par->SetLowerLimit(0);
			val = fLine->ampl;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/100);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Tail");
			tmp_par->SetLowerLimit(0);
			val = fLine->tail*fLine->ampl;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/50);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Sigma");
			tmp_par->SetLowerLimit(0);
			val = fLine->sigma;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/100);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Beta");
			tmp_par->SetLowerLimit(0);
			val = fLine->beta;
			tmp_par->SetValue(val);//Standard parametrization -> Using function peakFitFuncA(..)
			//tmp_par->SetValue(1/val);//Inverted parametrization -> Using function peakFitFuncB(..)
			tmp_par->SetMCMCStepSize(val/5);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Step");
			tmp_par->SetLowerLimit(0);
			val = fLine->step;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/10);
			//tmp_par->SetMCMCStepSize(fLine->ampl/1000);
			fParams.push_back(tmp_par);
			
			tmp_par = new Analysis::Param("Const");
			tmp_par->SetLowerLimit(0);
			val = fLine->cost;
			tmp_par->SetValue(val);
			tmp_par->SetMCMCStepSize(val/10);
			fParams.push_back(tmp_par);
			
			fParsN = fParams.size();
			
			fMaxPar = vector<double>(fParsN);
			if(fVerbosity>=1){
				cout << "-> " << fParsN << " parameters defined and initialised." << endl;
			}
			
			
			cout <<"\nMean init value: " << fLine->mean << endl;
			cout <<"Ampl init value: " << fLine->ampl << endl;
			cout <<"Tail init value: " << fLine->tail << endl;
			cout <<"Sigma init value: " << fLine->sigma << endl;
			cout <<"Beta init value: " << fLine->beta << endl;
			cout <<"Step init value: " << fLine->step << endl;
			cout <<"Constant init value: " << fLine->cost << endl;
			
			
			//Set the parameters step size for the MCMC algorithm
			if(!fParSteps){
				fParSteps = new vector<double>(fParsN);
			}else{
				if(fParSteps->size() != fParsN){
					delete fParSteps;
					fParSteps = new vector<double>(fParsN);
				}
			}
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				fParSteps->at(iPar) = fParams.at(iPar)->GetMCMCStepSize();
			}
			
			fExtParStep=true;
			
			fParInit = true;
			
			cout << "\nExiting from GammaLineLikelihood::DefineParameters()\n" << endl;
			
			return(0);
		}
		
		/*
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
		*/
		
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
			//Calculate the log-likelihood
			if(fEngineVerb>=4){
				cout << "\n====> Entering in GammaLineLikelihood::LogProb(...)." << endl;
			}
			double logprob =0;
			
			double xMin, xMax;
			
			int nBins = fHisto->GetNbinsX();
	
			for(int iBin=1; iBin<=nBins; iBin++){
				xMin = fHisto->GetBinLowEdge(iBin);
				xMax = fHisto->GetBinLowEdge(iBin+1);
				
				
				//Interpolate the function in the central interval and compare with the bin contents
				double lambda_i = ( peakFitFuncA(&xMax, (double*)(&(par.at(0))) ) + peakFitFuncA(&xMin, (double*)(&(par.at(0))) ) )*(xMax-xMin)/2.;
				
				//double lambda_i = peakFitFuncA( &xMid, (double*)(&(par.at(0))) );
				//double lambda_i = peakFitFuncB( &xMin, (double*)(&(par.at(0))) );
				
				double N_i = fHisto->GetBinContent(iBin);
				
				if(lambda_i>0){
					logprob += N_i*log(lambda_i) - lambda_i;
				}else{
					if( (((int)(N_i+0.5))>0) ) return -std::numeric_limits<double>::max();
				}
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
		
		
		double GammaLineLikelihood::MinLogProb(const double* parVal)
		{
			//This return the opposite of the LogProb used in the Minuit2 minimizer
	
			Analysis::GatorCalib::GammaLineLikelihood* this_ptr = static_cast<Analysis::GatorCalib::GammaLineLikelihood*>(Analysis::LikelihoodClassHolder::instance());
			
			vector<double> par(parVal, parVal + (this_ptr->fParsN) );
			
			return -( this_ptr->LogProb(par) );;
		}
		
		
		void GammaLineLikelihood::MetropolisMLE()
		{
			if(!fInit) return;
			
			Analysis::GatorCalib::GammaLineLikelihood::MetropolisMLE(NULL);
			fMinimized = true;
			return;
		}
		
		
		void GammaLineLikelihood::MetropolisMLE(vector<double> *initpar)
		{
			if(!fInit) return;
			
			cout << "Entering in LikelihoodClass::MetropolisMLE() function" << endl;
	
			//HERE starts the Metropilis search of MLE 
			vector<double> par(fParsN);
			vector<double> par_new(fParsN);
	
			vector<bool> fixedpars(fParsN);
			
			
			if(initpar){
				//Analysis::LikelihoodClass::SetInitPars(initpar);
				if(fInitPars) delete fInitPars;
				fInitPars = new vector<double>(fParsN);
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					fInitPars->at(iPar) = initpar->at(iPar);
				}
				fParInit = true;
				par = (*fInitPars);
			}
	
			//Initialize the parameters (to something reasonable)
			if(!fInitPars){
				if(fEngineVerb>=1){
					cout << "=> Self initialization of the parameters." << endl;
				}
				par =SelfParInit();
			}
			
			
			
			if(fEngineVerb>=1){
				cout << "=> Checking which are the fixed parameters." << endl;
			}
			unsigned nFixed = 0;
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				if(fParams.at(iPar)->IsFixed()){
					fixedpars.at(iPar) = true;
					par.at(iPar) = fParams.at(iPar)->GetValue();
					par_new.at(iPar) = fParams.at(iPar)->GetValue();
					nFixed += 1;
				}else{
					fixedpars.at(iPar) = false;
				}
			}
			if(fEngineVerb>=1){
				cout << "=> " << nFixed << " fixed parameters over " << fParsN << " total parameters." << endl;
			}
	
			if(nFixed>=fParsN){
				fMaxLogProb = LogProb(par);
				fMaxPar = par;
				if(fEngineVerb>=2){
					cout << "==> All Parameter fixed." << endl;
					cout << "==> Max LogProb: " << fMaxLogProb << endl;
					for(unsigned iPar=0; iPar<fParsN; iPar++){
						cout << "==>Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
					}
				}
				return;
			}
	
	
			double logprob = LogProb(par);
			double logprob_new;
	
			
			fMaxLogProb = logprob;
			fMaxPar = par;
	
			fRndGen->SetSeed();
			
			if(fEngineVerb>=1){
				cout << "=> Starting the Metropolis loop for the MLE finding with " << fNsteps << " steps and " << fMaxIter << " max iterations." << endl;
				
			}
			
			//Used to understand whther to increase or decrease the step size. Check every 500 trials.
			unsigned iStep=0, accepTrials=0, testTrials=0, iter=0;
			
			do{
				if(/*iStep%1000==0 &&*/ fEngineVerb>=2){
					cout << "==>Step: "<< iStep << endl;
					if(fEngineVerb>=3){
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							cout << "===> Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
						}
					}
					cout << "===> Logprob: " << logprob << endl;
				}
		
				bool accepted = false;
				MetrProp(par,par_new);
				
				
				if(/*iStep%1000==0 &&*/ fEngineVerb>=3){
					for(unsigned iPar=0; iPar<fParsN; iPar++){
						cout << "===> Proposed parameter <" << fParams.at(iPar)->GetName() << ">: " << par_new.at(iPar) << endl;
					}
				}
				
				//Check if the parameter is inside the allowed region if it is not in the range
				if(fEngineVerb>=3){
					cout << "===> Checking if new parameters are in range." << endl;
				}
				
				bool inRange=true;
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					if( fParams.at(iPar)->IsOutOfRange(par_new.at(iPar)) != 0 )
						inRange = false;
				}
				
				
				if(inRange){
					logprob_new = LogProb(par_new);
				
					if(/*iStep%1000==0 &&*/ fEngineVerb>=3){
						cout << "===> New LogProb:  " << logprob_new << endl;
					}
					
					if(logprob_new>=logprob){
						accepted = true;
					}else{
						double r = fRndGen->Rndm();
						if( (logprob_new-logprob)>log(r) ){
							accepted = true;
							if(fEngineVerb>=3){
								cout << "===> Accepted with odds = " << r << endl;
							}
						}else{
							if(fEngineVerb>=3){
								cout << "===> Rejected with odds = " << r << endl;
							}
						}
					}
				}else{
					if(/*iStep%1000==0 &&*/ fEngineVerb>=3){
						cout << "===> Rejected (out of range)" << endl;
					}
				}
				
				if( fEngineVerb>1 ){
					if(iter==0){
						cout << "==>Iteration 0 parameters status: " << endl;
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							cout << "=>  <" << fParams.at(iPar)->GetName() << "> old: " << par.at(iPar) << endl;
						}
						cout << "==> LogProb old: " << LogProb(par) << endl;
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							cout << "==>  <" << fParams.at(iPar)->GetName() << "> prop: " << par_new.at(iPar) << endl;
						}
						cout << "==> LogProb prop: " << LogProb(par_new) << endl;
						if(accepted){
							cout << "==>  Accepted" << endl;
						}else{
							cout << "==>  Rejected" << endl;
						}
					}
				}
				
				iter++;
				testTrials++;
				
				//Update step
				if(accepted){
					par = par_new;
					logprob = logprob_new;
			
					if(logprob>fMaxLogProb){
						fMaxLogProb = logprob;
						fMaxPar = par;
					}
			
					if(fTrackVals){
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							if(fParForced){
								fParams.at(iPar)->ForceValue(par.at(iPar),logprob);
							}else{
								fParams.at(iPar)->SetValue(par.at(iPar),logprob);
							}
						}
					}
					iStep++;
					accepTrials++;
				}
				
				if( (testTrials%500) == 0 ){
					//Check if the LogValue is a NAN and if the accepted trials are 0
					if( (accepTrials==0) && (LogProb(par)!=LogProb(par)) ){//Definition of a nan number
						//Increase the step size by factor 100 in order to go out of the bad region whith higher probability
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							fParSteps->at(iPar) = 100.*fParSteps->at(iPar);
						}
					}else{
						if( ((double)accepTrials) < (0.01*testTrials) ){
							if(fVerbosity>=1){
								cout << "-> Step: " << iStep << ", Iteration: " << iter << endl;
								if(testTrials>0) cout << "->  Acceptance frequency: " << (100.*accepTrials/testTrials) << "%. Reducing by factor 2 the size of the MCMC steps." << endl;
							
								if(fEngineVerb>=1){
									for(unsigned iPar=0; iPar<fParsN; iPar++){
										cout << "=>  Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
									}
									cout << "=>  LogProb: " << LogProb(par) << endl;
									for(unsigned iPar=0; iPar<fParsN; iPar++){
										cout << "=>  Proposed parameter <" << fParams.at(iPar)->GetName() << ">: " << par_new.at(iPar) << endl;
									}
									cout << "=>  Proposed LogProb: " << LogProb(par_new) << endl;
								}
								
							}
							for(unsigned iPar=0; iPar<fParsN; iPar++){
								fParSteps->at(iPar) = fParSteps->at(iPar)/2;
							}
						}
						
						if( ((double)accepTrials) > (0.15*testTrials) ){
							if(fVerbosity>=1){
								cout << "-> Step: " << iStep << ", Iteration: " << iter << endl;
								if(testTrials>0) cout << "->  Acceptance frequency: " << (100.*accepTrials/testTrials) << "%. Increasing by factor 2 the size of the MCMC steps." << endl;
							
								if(fEngineVerb>=1){
									for(unsigned iPar=0; iPar<fParsN; iPar++){
										cout << "=>  Parameter <" << fParams.at(iPar)->GetName() << ">: " << par.at(iPar) << endl;
									}
									cout << "=>  LogProb: " << LogProb(par) << endl;
									for(unsigned iPar=0; iPar<fParsN; iPar++){
										cout << "=>  Proposed parameter <" << fParams.at(iPar)->GetName() << ">: " << par_new.at(iPar) << endl;
									}
									cout << "=>  Proposed LogProb: " << LogProb(par_new) << endl;
								}
							
							}
							for(unsigned iPar=0; iPar<fParsN; iPar++){
								fParSteps->at(iPar) = fParSteps->at(iPar)*2;
							}
						}
					}
					
					
					testTrials = 0;
					accepTrials = 0;
				}
				
			}while( (iStep<fNsteps) && (iter<fMaxIter) );
			
			if(fEngineVerb>=1){
				cout << "=> End of Metropolis loop for the MLE finding." << endl;
				cout << "\n=> Max LogProb after MCMC: " << fMaxLogProb << endl;
			}
	
			cout << "Exiting from GammaLineLikelihood::MetropolisMLE() function." << endl;
	
			return;
		}
		
		
		void GammaLineLikelihood::Minuit2MLE(vector<double> &pars)
		{
			cout << "\nEntering in the GammaLineLikelihood::Minuit2MLE() function" << endl;
			
			
			if(fMinimizer){
				delete fMinimizer;
				fMinimizer=NULL;
			}
			if(fMinFunctor){
				delete fMinFunctor;
				fMinFunctor = NULL;
			}
			
			if(fVerbosity>1){
				cout << "-->Start the Minuit2 minimization with initial parameters:" << endl;
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					cout << "--> Parameter <" << fParams.at(iPar)->GetName() << ">: " << pars.at(iPar) << endl;
				}
			}
			
			
			//Register this class in the static class holder
			LikelihoodClassHolder::instance(this);
			//Now the MinLogProb static method will refer to this instance when using the class members
			
			
			fMinimizer = new ROOT::Minuit2::Minuit2Minimizer();
			
			// set tolerance , etc...
			fMinimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
			//fMinimizer->SetTolerance(0.001);
			//fMinimizer->SetPrecision(1e-10);
			fMinimizer->SetStrategy(2);
			
			//Setverbose in testing phase
			if(fVerbosity<=1){
				fMinimizer -> SetPrintLevel(-1);
			}
			if(fVerbosity>1){
				fMinimizer -> SetPrintLevel(0);
			}
			if(fVerbosity>2){
				fMinimizer -> SetPrintLevel(1);
			}
			
			
			//Set the instance before using the the static method MinLogProb in the minimizer
			Analysis::LikelihoodClassHolder::instance(this);
			//fMinFunctor = new ROOT::Math::Functor(this,&XurichAnalysis::XuLCElikelihood::MinLogProb,1);
			fMinFunctor = new ROOT::Math::Functor(&Analysis::GatorCalib::GammaLineLikelihood::MinLogProb, fParsN);
			
			fMinimizer->SetFunction(*fMinFunctor);
			
			
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				fMinimizer->SetLowerLimitedVariable(iPar, fParams.at(iPar)->GetName(), pars.at(iPar), fParSteps->at(iPar), 0.);
			}
			
			bool minFail = false;
			
			if(!fMinimizer->Minimize()){
				cout << "WARNING: Minimization in Minuit2 failed with error code " << fMinimizer->Status() << endl;
				minFail = true;
			}else{
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					fMaxPar.at(iPar) = (fMinimizer->X())[iPar];
				}
				
				//Minimize again without constrained parameters
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					fMinimizer->SetVariable(iPar, fParams.at(iPar)->GetName(), fMaxPar.at(iPar), fParSteps->at(iPar));
				}
				
				if(fEngineVerb>=2) cout << "=>Start the Minuit2 hessian calculation." << endl;
				
				if(!fMinimizer->Hesse()){
					cout << "WARNING: Hessian matrix in Minuit2 failed with error code " << fMinimizer->Status() << endl;
					return;
				}
			}
			
			if(minFail){
				string ans("");
				while(minFail){
					cout << "Do you want to run again the MCMC algoritm with the following initial parameters?" << endl;
					for(unsigned iPar=0; iPar<fParsN; iPar++){
						cout << " Parameter <" << fParams.at(iPar)->GetName() << ">: " << (fMinimizer->X())[iPar] << endl;
					}
					while(true){
						cout << "[Yes/No/Ignore]--> ";
						getline(cin,ans);
						if( (ans==string("y"))  || (ans==string("Y")) || (ans==string("n")) || (ans==string("N")) || (ans==string("I")) ) break;
					}
					
					if(ans==string("y") || ans==string("Y")){
						//Reset the parameters for a new run.
						vector<double> tmp_pars(fParsN);
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							tmp_pars.at(iPar) = (fMinimizer->X())[iPar];
							//Changing the parameters in order to find a better minimum
							fParSteps->at(iPar) = fParSteps->at(iPar)*10.;
						}
						fLogProbScaling = LogProb(tmp_pars);
						GammaLineLikelihood::MetropolisMLE(&tmp_pars);
						
						//Minimization using constrained parameters
						for(unsigned iPar=0; iPar<fParsN; iPar++){
							fMinimizer->SetLowerLimitedVariable(iPar, fParams.at(iPar)->GetName(), fMaxPar.at(iPar), fParSteps->at(iPar), 0.);
						}
						if(fMinimizer->Minimize()){
							minFail = false;
							for(unsigned iPar=0; iPar<fParsN; iPar++){
								tmp_pars.at(iPar) = (fMinimizer->X())[iPar];
							}
							fLogProbScaling = LogProb(tmp_pars);
							//Remove the lower bpundary for a more reliable fitting
							for(unsigned iPar=0; iPar<fParsN; iPar++){
								fMinimizer->SetVariable(iPar, fParams.at(iPar)->GetName(), tmp_pars.at(iPar), fParSteps->at(iPar));
							}
							fMinimizer->Minimize();
							
						}else{
							minFail = true;
							cout << "WARNING: Minimization in Minuit2 failed with error code " << fMinimizer->Status() << endl;
						}
					
					}else if( ans==string("I") ){
						//Use these parameters for the plot and for the chi2 and p-value calculations.
						minFail = false; //This exits from the while loop like if the fit converged.
					}else{
						fMinimized = false;
						cout << "Exiting from GammaLineLikelihood::Minuit2MLE(...) function." << endl;
						return;
					}
				}
			}
			
			
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				fMaxPar.at(iPar) = (fMinimizer->X())[iPar];
			}
			
			
			fLine->mean = (fMinimizer->X())[0];
			fLine->mean_err = sqrt(fMinimizer->CovMatrix(0,0));
			
			fLine->ampl = (fMinimizer->X())[1];
			fLine->ampl_err = sqrt(fMinimizer->CovMatrix(1,1));
			
			fLine->tail = (fMinimizer->X())[2];
			fLine->tail_err = sqrt(fMinimizer->CovMatrix(2,2));
			
			fLine->sigma = (fMinimizer->X())[3];
			fLine->sigma_err = sqrt(fMinimizer->CovMatrix(3,3));
			
			fLine->beta = (fMinimizer->X())[4];
			fLine->beta_err = sqrt(fMinimizer->CovMatrix(4,4));
			
			fLine->step = (fMinimizer->X())[5];
			fLine->step_err = sqrt(fMinimizer->CovMatrix(5,5));
			
			fLine->cost = (fMinimizer->X())[6];
			fLine->cost_err = sqrt(fMinimizer->CovMatrix(6,6));
			
			fLine->chi2 = GetChi2();
			fLine->chi2ndof = GetChi2NDof();
			
			fLine->p_value = GetPval();
			fLine->p_value_ndof = GetPvalNDof();
			
			
			
			cout << "\nMinuit2 MLE results: "  << endl;
			cout << "Min LL value = " << fMinimizer->MinValue() << endl;
			for(unsigned iPar=0; iPar<fParsN; iPar++){
				cout << fParams.at(iPar)->GetName() << " MLE = (" << (fMinimizer->X())[iPar] << " +- " << sqrt(fMinimizer->CovMatrix(iPar,iPar)) << ")" << endl;
			}
			cout << "Chi2 = " << fLine->chi2 << endl;
			cout << "Chi2/ndf = " << fLine->chi2ndof << endl;
			cout << "p-value = " << fLine->p_value_ndof << endl;
			
			
			
			cout << "Exiting from GammaLineLikelihood::Minuit2MLE(...) function." << endl;
			
			fMinimized = true;
			
			
			return;
		}
		
		
		void GammaLineLikelihood::Minuit2MLE()
		{
			if(!fMinimized) return;
			
			GammaLineLikelihood::Minuit2MLE(fMaxPar);
			return;
		}
		
		
		void GammaLineLikelihood::CalculatePValueLikelihood(vector<double> *pars)
		{
			if(fVerbosity<=1){
				cout << "GammaLineLikelihood::CalculatePValueLikelihood() function" << endl;
			}
			
			double *parsArr = new double[fParsN];
			
			if(pars){
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					parsArr[iPar] = pars->at(iPar);
				}
			}else{
				if(!fMinimized) return;
				for(unsigned iPar=0; iPar<fParsN; iPar++){
					parsArr[iPar] = fMaxPar.at(iPar);
				}
			}
			
			
			
		// initialize test statistic -2*lambda
			double logLambda = 0.0;
			double xMin, xMax, lambda_i;
			int nBins = fHisto->GetNbinsX();
			
			for(int iBin = 1; iBin <= fNbins; ++iBin) {
				// get the number of observed events
				double n_i = fHisto->GetBinContent(iBin);
				
				// get the number of expected events using the interpolation of the function in the bin "iBin"
				xMin = fHisto->GetBinLowEdge(iBin);
				xMax = fHisto->GetBinLowEdge(iBin+1);
				
				double lambda_i = (peakFitFuncA(&xMax, parsArr ) + peakFitFuncA(&xMin, parsArr))*(xMax-xMin)/2.;
				//lambda_i = peakFitFuncA( &xMin, parsArr );
				//lambda_i = peakFitFuncB( &xMin, parsArr );
				
				// get the contribution from this datapoint
				if(lambda_i <= 0.0){
					if( ((int)(n_i+0.5))>0 ){
						delete [] parsArr;
						return;
					}
				}
				if (n_i == 0){
					logLambda += lambda_i;
				}else{
					logLambda += (lambda_i - n_i) + n_i * log(n_i / lambda_i);
				}
			}

			// rescale
			logLambda *= 2.0;

			//p value from chi^2 distribution, returns zero if logLambda < 0
			fPval = TMath::Prob( logLambda, fNbins );
			fPvalNDof = TMath::Prob( logLambda, GetNDoF() );
			
			fPvalues = true;
			
			delete [] parsArr;
			
			if(fVerbosity<=1){
				cout << "Exiting from GammaLineLikelihood::Minuit2MLE(...) function." << endl;
			}
			
			// no error
			return;
		}
		
		
		double GammaLineLikelihood::GetPval(vector<double> *pars)
		{
			if(pars){
				if(pars->size()!=fParsN){
					cout << "WARNING --> GammaLineLikelihood::GetPval(...): the vector \"pars\" has wrong size! Returning -1." << endl;
					return -1.0;
				}
				CalculatePValueLikelihood(pars);
				return fPval;
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if(!fPvalues){
				cout << "WARNING --> GammaLineLikelihood::GetPval(...): the calculation of the pvalues failed! Returning -1." << endl;
				return -1.0;
			}
			
			return fPval;
		}
		
		
		double GammaLineLikelihood::GetPvalNDof(vector<double> *pars)
		{
			if(pars){
				if(pars->size()!=fParsN) return -1.0;
				CalculatePValueLikelihood(pars);
				return fPvalNDof;
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if(!fPvalues) return -1.0;
			
			return fPvalNDof;
		}
		
		
		double GammaLineLikelihood::GetChi2(vector<double> *pars)
		{
			if(pars){
				if(pars->size()!=fParsN){
					cout << "WARNING --> GammaLineLikelihood::GetChi2(...): the vector \"pars\" has wrong size! Returning -1." << endl;
					return -1.0;
				}
				CalculatePValueLikelihood(pars);
				if( !( (fPvalNDof>0.0) && (fPvalNDof<1.0) ) ){
					cout << "WARNING --> GammaLineLikelihood::GetChi2(...): the calculation of the pvalues failed! Returning -1." << endl;
					return -1.0;
				}
				return ROOT::Math::chisquared_quantile_c( fPvalNDof, (double)GetNDoF() );
			}
			
			if(!fPvalues) CalculatePValueLikelihood();
			
			if( !( (fPvalNDof>0.0) && (fPvalNDof<1.0) ) ){
				cout << "WARNING --> GammaLineLikelihood::GetChi2(...): the calculation of the pvalues failed! Returning -1." << endl;
				return -1.0;
			}
			
			return ROOT::Math::chisquared_quantile_c( fPvalNDof, (double)GetNDoF() );
		}
		
		
		double GammaLineLikelihood::GetChi2NDof(vector<double> *pars)
		{
			if(pars){
				if(pars->size()!=fParsN) return -1.0;
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



