#ifndef GATOR_STRUCTS_H
#define GATOR_STRUCTS_H

#include <string>

using namespace std;

typedef struct LineStruct{
	
	string linename; //Mass+element+energy
	
	string mass; //Put here also whther it is a metastable isomer
	string element;
	double Energy;
	double Energy_err;
	double BR;
	double BR_err;
	double Effic;
	double Effic_err;
	double BRxEffic;
	double BRxEffic_err;
	
	//Variables used to define the ideal limits of the three counting regions (the units can be either energy or MCA channels)
	//Sample
	double MCAcenter; //From the calibration
	double MCAsigma; //From the calibration
	double LowEdgeL; //Low edge of the left compton region
	double UpEdgeL; //Up edge of the left compton region
	double LowEdgeS; //Low edge of the signal region
	double UpEdgeS; //Up edge of the signal region
	double LowEdgeR; //Low edge of the right compton region
	double UpEdgeR; //Up edge of the right compton region
	//Background
	double BgMCAcenter; //From the calibration
	double BgMCAsigma; //From the calibration
	double BgLowEdgeL; //Low edge of the left compton region
	double BgUpEdgeL; //Up edge of the left compton region
	double BgLowEdgeS; //Low edge of the signal region
	double BgUpEdgeS; //Up edge of the signal region
	double BgLowEdgeR; //Low edge of the right compton region
	double BgUpEdgeR; //Up edge of the right compton region
	
	
	//This variables are for the actual fits
	
	//Signal
	double binL; //Number of bins at the left of the signal zone
	double countsL; //Number of counts found in the left zone
	double binS; //Number of bins in the central signal zone
	double countsS; //Number of counts found in the signal zone
	double binR; //Number of bins at the right of the signal zone
	double countsR; //Number of counts found in the right zone
	
	//Background
	double BgbinL; //Number of bins at the left of the signal zone
	double BgcountsL; //Number of counts found in the left zone
	double BgbinS; //Number of bins in the central signal zone
	double BgcountsS; //Number of counts found in the signal zone
	double BgbinR; //Number of bins at the right of the signal zone
	double BgcountsR; //Number of counts found in the right zone
	
	
	//Variables used only by the "sampleanalysis" program
	//The line
	double RunTime;
	double PeakCnts; //This should represent the estimated total (gross) rate in the signal region
	double PeakCnts_err;
	double CompCnts; //Here there should be the estimated compton counts in the signal region
	double CompCnts_err;
	double LineCnts; //This should represent the estimated signal (net) rate
	double LineCnts_err;
	double LdCnts; //Detection limit for the activity (in counts)
	bool Detected;
	
	//The background
	double BgRunTime;
	double BgPeakCnts; //This should represent the estimated total (gross) rate in the signal region
	double BgPeakCnts_err;
	double BgCompCnts; //Here there should be the estimated compton counts in the signal region
	double BgCompCnts_err;
	double BgLineCnts; //This should represent the estimated signal (net) rate
	double BgLineCnts_err;
	double BgLdCnts; //Detection limit for the background in counts
	bool BgDetected;
	
	//Final activity
	double LdActiv; //Detection limit for the activity (in mBq/kg)
	double Activ; //Activity calculated anyway (if there is no detection this number doesn't have any meaning)
	double Activ_err; //Activity calculated anyway (if there is no detection this number doesn't have any meaning)
}LineStruct;


typedef struct ResultsStruct{
	double RunTime;
	double EffActiv;
	double EffActiv_err;
	double RelActivErr;
	int DetecNum;
	double EffDetLim;
	double ScaleFact;
}ResultsStruct;


typedef struct CalibLine{
	
	string massN;
	string element;
	double litEn;
	double litEn_err;
	double MCAlowch;
	double MCAupch;
	
	double mean;
	double mean_err;
	double sigma;
	double sigma_err;
	double beta;
	double beta_err;
	double ampl;
	double ampl_err;
	double tail;
	double tail_err;
	double ratio;
	double ratio_err;
	
	double step;
	double step_err;
	double cost;
	double cost_err;
	
	double p_value;
	double p_value_ndof;
	
}CalibLine;


typedef struct RateStruct{
	
	double totcounts;
	double runtime;//This can be days or seconds..... it is assigned in the script
	double starttime;//This can be days or seconds..... it is assigned in the script
	
}RateStruct;


#endif
