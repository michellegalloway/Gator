#include "GatorParam.hh"

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


using namespace std;


GatorParam::GatorParam(string name, double upEdge, bool trackvals):
    fName(name),
	fLowerLimit(0.),
	fUpperLimit(upEdge),
	fHist(NULL),
	fFixed(false),
	fValue(0),
	fLogProb(0),
	fNbins(1000),
	fHistoFlag(false),
	fTrackVals(trackvals)
{
	
	fLockLowEdge = true;
	//Create 2 empty vectors (the likelihood value and )
	fValsChain = vector<vector<double> >(2);
}


GatorParam::GatorParam(string name, double lowEdge, double upEdge, bool trackvals):
    fName(name),
	fLowerLimit(lowEdge),
	fUpperLimit(upEdge),
	fHist(NULL),
	fFixed(false),
	fValue(0),
	fLogProb(0),
	fNbins(1000),
	fHistoFlag(false),
	fTrackVals(trackvals)
{
	
	fLockLowEdge = false;
	//Create 2 empty vectors (the parameter value and the likelihood value)
	fValsChain = vector<vector<double> >(2);
	
	if(lowEdge<0) fLockLowEdge=false;
	
}


void GatorParam::SetLimits(double lowerlimit, double upperlimit){
	if(lowerlimit>upperlimit){
		cerr << '\nGatorParam::SetLimits(...) --> ERROR: Tryng to set lower edge higher than upper edge for parameter "' << fName << '"' << endl << endl;
		return;
	}
	
	if(fLockLowEdge){
		if(lowerlimit<0){
			fLowerLimit=0;
		}else{
			fLowerLimit=lowerlimit;
		}
	}else{
		fLowerLimit=lowerlimit;
	}
	
	if(upperlimit>fLowerLimit){
		fUpperLimit = upperlimit;
	}else{
		fUpperLimit=fLowerLimit;
		fFixed = true;
	}
	
	
}


bool GatorParam::SetValue(double value, double logprob){
	
	bool flag=true;
	if(value < fLowerLimit){
		fValue = fLowerLimit;
		flag = false;
		//return false;
	}
	if(value > fUpperLimit){
		fValue = fUpperLimit;
		flag = false;
		//return false;
	}
	if(flag){
		fValue = value;
		fLogProb = logprob;
	}
	
	
	if(fTrackVals){
		fValsChain.at(0).push_back(fValue);
		fValsChain.at(1).push_back(fLogProb);
	}
	
	return flag;
	
}


void GatorParam::ForceValue(double value, double logprob){
	
	
	fValue = value;
	fLogProb = logprob;
	
	
	if(fTrackVals){
		fValsChain.at(0).push_back(fValue);
		fValsChain.at(1).push_back(fLogProb);
	}
	
	return;
	
}


int GatorParam::IsOutOfRange(double value) const{
	
	int retval;
	
	if( (fLowerLimit <= value) && (value <= fUpperLimit) ) retval =0;
	
	if( (fLowerLimit > value) ) retval =-1;
	if( (fUpperLimit < value) ) retval =1;
	
	return retval;
	
}


void GatorParam::Fix(double value){
	if(fTrackVals&&(fValsChain.size()>0)){
		//cout << "\nGatorParam::Fix(...) cannot fix the parameter:  " << fName << endl;
		fTrackVals = false; //When I fix the parameter I don't track it any more!
	}
	
	fFixed = true;
	fValue = value;
	
	//cout << '\nFixing "'  << fName << '" parameter to ' << fValue;
	return;
}