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


#include "GatorLikelihoodClass.hh"
#include "GatorParam.hh"



using namespace std;


GatorLikelihoodClass::GatorLikelihoodClass( string sourcefilename, string treename ):
	GatorTreeDataLoader(sourcefilename,treename),
	fDataType(kTreeData),
	m_nPar(0),
	fLinesInit(false),
	fParInit(false),
	fMaxLogProb(0.),
	fWidthScale(3.) //That's a Gator standard for the analisys
{
	
	mb_BGflag = false;
	
	//gROOT->ProcessLine("#include <vector>");
	
	fRndGen = new TRandom3(0);
	
	SetPoissonCounting();
	
	fMinuit = NULL;
	
	//Initialization of the arglist used for the MIGRAD routine.
	//To have a better understanding make an interactive session in root and use TMinuit::mnhelp("MIGrad")
	//The first argument for MIGRAD is the maximum number of function calls
	fMinuitArglist[0] = 1000; //Optional
	//This is the tollerance
	fMinuitArglist[1] = 0.01; //Optional
	
	//w = NULL;
	Lcnts = NULL;
	Gcnts = NULL;
	Rcnts = NULL;
	ParSteps = NULL;
	//effic_sum = 0;
	
}





GatorParam* GatorLikelihoodClass::GetParameter(unsigned iPar){
	if(iPar<2*m_nPar+1){
		return fParams.at(iPar);
	}else{
		cerr << "\nOut of range: parameter index " << iPar << " doesn't exist. Abort program!\n" << endl;
		return NULL;
	}
};




//THIS ALGORITHM DOESN'T WORK........ PORCAMADONNA!!!
/*
void GatorLikelihoodClass::MLEMultiLine(){
	
	cout << "Entering in GatorLikelihoodClass::MLEMultiLine() function" << endl;
	
	if(!fInitialized){
		cerr << "\GatorLikelihoodClass::MLEMultiLine() --> ERROR: The class is not properly initialized! Returning a non sense value!" << endl;
		return;
	}
	
	double logprobcomp, logprobcomp_new;
	double logprob, logprob_new;
		
	vector<double> par(m_nPar);
	vector<double> par_new(m_nPar);
	
	vector<double> S_est(m_nLines);
	
	//Initialize the nuisance parameters (expected background)
	for(unsigned iLine=0; iLine<m_nLines; iLine++){
		par.at(2*iLine+1) = Lcnts->at(iLine);
		par.at(2*iLine+2) = Rcnts->at(iLine);
	}
	
	double diff;
	double fact;
	
	
	
	if( !(fParams.at(0)->IsFixed()) ){//Here the activity is varied (finding the global MLE)
		
		//Initialize the activity
		par.at(0) = 0;
		
		double effic_sum = 0; //I will need it later for the Estimation step
		
		for(unsigned iLine=0; iLine<m_nLines; iLine++){
			
			par.at(0) += Gcnts->at(iLine) - w->at(iLine)*(Lcnts->at(iLine)+Rcnts->at(iLine))/(m_nLines*(mv_lines->at(iLine)).BRxEffic);
			
		}
		
		logprobcomp = LogProbComplete(par); 
		logprob = LogProb(par);
		
		unsigned step = 0;
		
		//Start the Expectation-Maximization loop
		do{
			cout << "\nStep:  "<< step << "    Signal:  " << par.at(0) << "    Logprob:  " << logprobcomp << endl;
			
			par_new.at(0) = 0; //The new activity is recomputed at each loop
			
			//Estimation step with the (i-1) values of parameters
			for(unsigned iLine=0; iLine<m_nLines; iLine++){
				
				S_est.at(iLine) = par.at(0)*(mv_lines->at(iLine)).BRxEffic/( par.at(0)*(mv_lines->at(iLine)).BRxEffic + w->at(iLine)*(par.at(2*iLine+1)+par.at(2*iLine+2)) );
				
				//Here I am doing also the Maximization step
				par_new.at(0) += S_est.at(iLine)/effic_sum;
				
				fact = (Gcnts->at(iLine)+Lcnts->at(iLine)+Rcnts->at(iLine)-S_est.at(iLine))/(1+w->at(iLine))/(Lcnts->at(iLine)+Rcnts->at(iLine));
				
				par_new.at(2*iLine+1) = fact*(Lcnts->at(iLine));
				par_new.at(2*iLine+2) = fact*(Rcnts->at(iLine));
				
			}
			
			
			logprobcomp_new = LogProbComplete(par_new);
			logprob_new = LogProb(par_new);
			
			bool ParIsChanged = false;
			for(unsigned iPar=0; iPar<m_nPar; iPar++){
				int ecode = fParams.at(iPar)->IsOutOfRange(par_new.at(iPar));
				if(ecode==-1){
					par_new.at(iPar) = fParams.at(iPar)->GetLowerLimit();
					ParIsChanged = true;
				}
				if(ecode==+1){
					par_new.at(iPar) = fParams.at(iPar)->GetUpperLimit();
					ParIsChanged = true;
				}
			}
			
			//Recompute the complete logprob if one of the parameters has changed after the optimization
			if(ParIsChanged) {
				logprobcomp_new = LogProbComplete(par_new);
				logprob_new = LogProb(par_new);
			}
			
			//diff = logprobcomp_new-logprobcomp;
			diff = logprob_new-logprob;
			
			if(diff<0){
				cout << '\GatorLikelihoodClass::MLEMultiLine() --> WARNING: The updated "completed likelihood" at iteration' << step << ' is lower than the old value.' << endl << endl;
			}
			
			for(unsigned iPar=0; iPar<m_nPar; iPar++){
				fParams.at(iPar)->SetValue(par_new.at(iPar),logprob);
				//fParams.at(iPar)->SetValueForce( par.at(iPar),logprob );
				par.at(iPar) = par_new.at(iPar);
			}
			
			
			logprobcomp = logprobcomp_new;
			logprob = logprob_new;
			
			if( TMath::Abs(diff)<=fAbsPrec ){
				//cout << "Step " << step << "  -->  logprob abs diff: " << TMath::Abs(diff) << endl;
				break;
			}
			
			step++;
			
		}while(step<fMaxIter);
		
		
		
	}else{//Here I keep the activity fixed in order to compute the profile likelihood
		
		par.at(0) = fParams.at(0)->GetValue();
		
		for(unsigned iLine=0; iLine<m_nLines; iLine++){
			
			S_est.at(iLine) = par.at(0)*(mv_lines->at(iLine)).BRxEffic/( par.at(0)*(mv_lines->at(iLine)).BRxEffic + w->at(iLine)*(par.at(2*iLine+1)+par.at(2*iLine+2)) );
			
		}
		
		logprobcomp = LogProbComplete(par); 
		
		unsigned step = 0;
		
		par_new.at(0) = par.at(0);
		
		do{
			//cout << "\nStep:  "<< step << "    Signal:  " << par.at(0) << "    Logprob:  " << logprob << endl;
			
			//Estimation step with the (i-1) values of parameters
			for(unsigned iLine=0; iLine<m_nLines; iLine++){
				
				S_est.at(iLine) = par.at(0)*(mv_lines->at(iLine)).BRxEffic/( par.at(0)*(mv_lines->at(iLine)).BRxEffic + w->at(iLine)*(par.at(2*iLine+1)+par.at(2*iLine+2)) );
				
				//Here I am doing also the Maximization step
				fact = (Gcnts->at(iLine)+Lcnts->at(iLine)+Rcnts->at(iLine)-S_est.at(iLine))/(1+w->at(iLine))/(Lcnts->at(iLine)+Rcnts->at(iLine));
				
				par_new.at(2*iLine+1) = fact*(Lcnts->at(iLine));
				par_new.at(2*iLine+2) = fact*(Rcnts->at(iLine));
				
			}
			
			
			//Recalculate the likelihood
			logprobcomp_new = LogProbComplete(par_new);
			logprob_new = LogProb(par_new);
			
			bool ParIsChanged = false;
			for(unsigned iPar=0; iPar<m_nPar; iPar++){
				if(!fParams.at(iPar)->IsFixed()){
					int ecode = fParams.at(iPar)->IsOutOfRange(par_new.at(iPar));
					if(ecode==-1){
						par_new.at(iPar) = fParams.at(iPar)->GetLowerLimit();
						ParIsChanged = true;
					}
					if(ecode==+1){
						par_new.at(iPar) = fParams.at(iPar)->GetUpperLimit();
						ParIsChanged = true;
					}
				}else{
					par_new.at(iPar) = fParams.at(iPar)->GetValue();//Restore the fixed value
					ParIsChanged = true;
				}
			}
			
			
			if(ParIsChanged) {
				logprobcomp_new = LogProbComplete(par_new);
				logprob_new = LogProb(par_new);
			}
			
			//diff = logprobcomp_new-logprobcomp;
			diff = logprob_new-logprob;
			
			if(diff<0){
				cout << '\GatorLikelihoodClass::MLEMultiLine() --> WARNING: The updated "completed likelihood" at iteration' << step << ' is lower than the old value.' << endl << endl;
			}
			
			//Set the values for every parameter but the activity (par[0])
			for(unsigned iPar=1; iPar<m_nPar; iPar++){
				fParams.at(iPar)->SetValue(par_new.at(iPar),logprob);
				//fParams.at(iPar)->SetValueForce( par.at(iPar),logprob );
				par.at(iPar) = par_new.at(iPar);
			}
			
			logprobcomp = logprobcomp_new;
			logprob = logprob_new;
			
			if( TMath::Abs(diff)<=fAbsPrec ){
				//cout << "Step " << step << "  -->  logprob abs diff: " << TMath::Abs(diff) << endl;
				break;
			}
			
			step++;
			
		}while(step<fMaxIter);
		
	}
	
	cout << "Exiting from GatorLikelihoodClass::MLEMultiLine() function" << endl;
	
	return;
}

//THIS ALGORITHM DOESN'T WORK
void GatorLikelihoodClass::MLEMultiLineNoActiv(){
	
	fParams.at(0)->Fix(0);
	 MLEMultiLine();//THIS ALGORITHM DOESN'T WORK
	
	return;
}*/


/*
vector<double> GatorLikelihoodClass::GetProfileLikelihoods(double snr){
	
	cout << "\nEntering in GatorLikelihoodClass::GetProfileLikelihoods(...)" << endl;
	
	vector<double> likelihoodvec;
	
	double sign = snr*fSigma.at(0)*fBgLevel.at(0);
	
	cout << " snr: " << snr << " --> signal: " << sign << endl;
	
	fParams.at(0)->Fix(sign);
	
	for(unsigned iEntry=0; iEntry<fTreeEntries ;iEntry++ ){
		
		AnalyzeTreeEntry(iEntry);
		
		likelihoodvec.push_back(fMaxLogProb);
		
	}
	
	cout << "Exiting from GatorLikelihoodClass::GetProfileLikelihoods(...)\n" << endl;
	
	return likelihoodvec;
}*/
