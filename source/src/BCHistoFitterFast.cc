#include "BCHistoFitterFast.hh"

#include "BAT/BCDataSet.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCMath.h"


using namespace std;

Gator::BcHistoFitterFast::BcHistoFitterFast(TH1D * hist, TF1 * func):BCHistogramFitter(hist,func)
{
	if( !(hist && func) ) exit(-1);
	
	fFitFunction = func;
	fHistogram = hist;
	
	fFlagIntegration = false;
	
	nBins = fHistogram->GetNbinsX();
	fHistVec.resize(nBins);
	fLowEdges.resize(nBins);
	fUpEdges.resize(nBins);
	
	for(int iBin=1; iBin<=nBins; iBin++)
	{
		fHistVec.at(iBin-1) = fHistogram->GetBinContent(iBin);
		fLowEdges.at(iBin-1) = fHistogram->GetBinLowEdge(iBin);
		fUpEdges.at(iBin-1) = fHistogram->GetBinLowEdge(iBin+1);
	}
}


double Gator::BcHistoFitterFast::LogLikelihood(const std::vector<double> & parameters)
{
	// initialize probability
	double loglikelihood = 0;

	// set the parameters of the function
	fFitFunction->SetParameters(&parameters.at(0));
	
	// loop over all bins
	for (int iBin=0; iBin<nBins; ++iBin)
	{
		// get function value at lower bin edge
		double fedgelow = fFitFunction->Eval(fLowEdges.at(iBin));
		
		// get function value at upper bin edge
		double fedgehi = fFitFunction->Eval(fUpEdges.at(iBin));
		
		
		// get the number of observed events
		double y = fHistVec.at(iBin);

		double lambda = 0.;

		// use ROOT's TH1D::Integral method
		if(fFlagIntegration)
		{
			lambda = fFitFunction->Integral(fLowEdges.at(iBin), fUpEdges.at(iBin));
		}
		else
		{// use linear interpolation
			lambda = (fUpEdges.at(iBin)-fLowEdges.at(iBin)) * (fedgehi + fedgelow)/2;
		}

		// get the value of the Poisson distribution
		loglikelihood += BCMath::LogPoisson(y, lambda);
	}

	return loglikelihood;
}