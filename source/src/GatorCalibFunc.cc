#include "GammaLineLikelihood.hh"
#include "GatorStructs.h"
#include "GatorCalibFitters.h"
#include "screenfncs.h"
#include "misc.h"
#include "BCHistoFitterFast.hh"

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>
#include <BAT/BCParameter.h>
#include <BAT/BCH2D.h>

#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "Math/Functor.h"

#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>



using namespace std;

vector<string>& loadConfFile(string filename)
{
	
	vector<string>* lines_vec = new vector<string>;
	
	if(!fexist(filename)){
		//ERROR MESSAGE
		return *lines_vec; 
	}
	
	string line;
	
	ifstream infile(filename.c_str());
	
	getline(infile,line); //Throw out the first line (it is an header)
	
	while(getline(infile,line)){
		lines_vec->push_back(line);
	}
	
	return *lines_vec;
}


CalibLine lineInit(string line)
{
	CalibLine calibline;
	
	stringstream linestream("");
	
	linestream << line;
	
	linestream >> line;
	calibline.massN = line;
	
	linestream >> line;
	calibline.element = line;
	
	linestream >> line;
	calibline.litEn = atof(line.c_str());
	
	linestream >> line;
	calibline.litEn_err = atof(line.c_str());
	
	linestream >> line;
	calibline.MCAlowch = atof(line.c_str());
	
	linestream >> line;
	calibline.mean = atof(line.c_str());
	
	linestream >> line;
	calibline.MCAupch = atof(line.c_str());
	
	linestream >> line;
	calibline.sigma = atof(line.c_str());
		
	linestream >> line;
	calibline.beta = atof(line.c_str()); 
	
	linestream >> line;
	calibline.tail = atof(line.c_str());
	
	calibline.histo = NULL;
	
	return calibline;
	
}


extern TApplication* theApp;

bool doFitBAT(TH1D* MCAhisto, CalibLine& line)
{
	//Make a copy of the ROI of the original histo and put it in a new histogram
	int firstbin = MCAhisto->FindBin(line.MCAlowch);
	int lastbin = MCAhisto->FindBin(line.MCAupch);
	int nbins = 1 + lastbin - firstbin;
	double xmin = MCAhisto->GetBinLowEdge(firstbin);
	double xmax = MCAhisto->GetBinLowEdge(lastbin+1);
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << line.element << line.massN << "_" << (int)(line.litEn+0.5) << "keV" ;
	TH1D* tmphisto = new TH1D(ss_histoname.str().c_str(),";Channel;Counts",nbins,xmin,xmax);
	//tmphisto -> SetDirectory(0);
	
	
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = MCAhisto->GetBinCenter(bin);
		tmphisto->Fill(bincenter,MCAhisto->GetBinContent(bin));
	}
	
	CalibLine intialValues = line;//Save the initial parameters as read from the file
	
	
	//Start with the fitting process
	TF1* ff_MCA = new TF1( "ff_MCA", peakFitFuncB, xmin, xmax, 7 );
	ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
	//ff_MCA->SetParNames("Mean","Ampl","Ratio","Sigma","Beta","Constant");
	
	
	int nPars = ff_MCA->GetNpar();
	
	vector<bool> fixedpars(nPars,false);
	
	bool tail_flag = true;
	bool reuse_flag = false;
	
	START:
	
	string ans;
	
	if(!reuse_flag){
		while(true){
			cout << "\nDo you want to use the low energy tail for the fit?\n" << "[y/n] > ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
		}

		if(ans==string("y") || ans==string("Y")){
			tail_flag = true;
		}else{
			tail_flag = false;
		}
	}else{
		tail_flag = true;
	}
	
	
	
	//This must be initialized with this order
	if(!reuse_flag){
		costInit(tmphisto,line);
		stepInit(tmphisto,line);
		amplInit(tmphisto,line);
		//sigmaInit(tmphisto,line);
	}
	
	
	
	
	//Set initial parameters
	ff_MCA -> SetParameter(0, line.mean);
	ff_MCA -> SetParameter(1, line.ampl);
	if(tail_flag){
		ff_MCA -> SetParameter(2, line.tail);
		fixedpars.at(2) = false;
	}else{
		ff_MCA -> FixParameter(2, 0.0);
		fixedpars.at(2) = true;
	}
	ff_MCA -> SetParameter(3, line.sigma);
	if(tail_flag){
		ff_MCA -> SetParameter(4, line.beta);
		fixedpars.at(4) = false;
	}else{
		ff_MCA -> FixParameter(4, 0.0);
		fixedpars.at(4) = true;
	}
	if(line.step>0.0){
		ff_MCA -> SetParameter(5, line.step);
	}else{
		ff_MCA -> SetParameter(5, 1e-80);
	}
	
	ff_MCA -> SetParameter(6, line.cost);
	
	
	
	
	
	if(reuse_flag){//It means that the previous fit iteration was performed without the tail
		double val, err, min, max;
		
		if( (line.mean_err/line.mean)<1e-3 ){
			//If the error is too small I allow for a variation of 0.1% around the the previous fit value
			ff_MCA -> SetParLimits(0, 0.999*line.mean, 1.001*line.mean );
		}else{
			min = TMath::Max(line.mean-3*line.mean_err, xmin);
			max = TMath::Min(line.mean+3*line.mean_err, xmax);
			ff_MCA -> SetParLimits(0, min, max ); //Mean
		}
		
		{//The amplitude is allowed to variate for a large range determined from the previous fit, but cannot get the 0 value!
			min = TMath::Max(line.ampl-3*line.ampl_err, 1e-80);
			ff_MCA -> SetParLimits(1, min, line.ampl+3*line.ampl_err ); //Gaussian amplitude
		}
		
		ff_MCA -> SetParLimits(2, 1e-80, 1.); //Tail ratio
		
		{//The width is allowed to variate for a 10% around it's previous fit value
			ff_MCA -> SetParLimits(3, 0.9*line.sigma, 1.1*line.sigma ); //Sigma
		}
		
		ff_MCA -> SetParLimits(4, 1e-80, 10*line.beta ); //Beta (or Gamma when parametrization B is used)
		
		if( !(line.step>0.0) ){
			min = TMath::Max(line.step-3*line.step_err, 1e-80);
			ff_MCA -> SetParLimits(5, min, line.step+3*line.step_err ); //Step amplitude	
		}else{
			ff_MCA -> SetParLimits(5, 1e-80, line.ampl/2 );//Step amplitude
		}
		
		if(line.cost_err/line.cost<1e-2){
			//The constant is allowed to variate for a 1% of the previous fit if the error is tto small
			ff_MCA -> SetParLimits(6, 0.99*line.cost, 1.01*line.cost ); //Constant
		}else{
			ff_MCA -> SetParLimits(6, line.cost-3*line.cost_err, line.cost+3*line.cost_err ); //Constant
		}
		
		
	}
	else{
		
		ff_MCA -> SetParLimits(0, 0.99*line.mean, 1.01*line.mean ); //Gaussian center
		
		ff_MCA -> SetParLimits(1, 1e-80, 10*line.ampl ); //Gaussian Amplitude
		
		//ff_MCA -> SetParLimits(2, 0.0, std::numeric_limits<double>::infinity() ); //Tail strenght
		
		if(tail_flag) ff_MCA -> SetParLimits(2, 1e-80, 1.); //Tail ratio
		
		ff_MCA -> SetParLimits(3, 0.9*line.sigma, 1.1*line.sigma );
		
		if(tail_flag) ff_MCA -> SetParLimits(4, 1e-80, 10*line.beta ); //Beta
		
		if( line.step>0.0 ){
			ff_MCA->SetParLimits(5, 1e-80, 10*line.step );//Step
		}else{
			ff_MCA->SetParLimits(5, 1e-80, line.ampl/2 );//Step
		}
		
		ff_MCA -> SetParLimits(6, 1e-80, 10*line.cost );
	}
	
	
	
	//ff_MCA -> FixParameter(0, 0.0 );
	//ff_MCA -> FixParameter(1, 0.0 ); //Gaussian ampl
	//ff_MCA -> FixParameter(2, 0.0 ); //Tail ampl
	//if(!tail_flag) ff_MCA -> FixParameter(2, 0.0 ); //Tail ratio
	//ff_MCA -> FixParameter(3, 0.0 );
	//if(!tail_flag) ff_MCA -> FixParameter(4, 0.0 ); //Beta
	//ff_MCA -> FixParameter(5, 0.0 );
	//ff_MCA -> FixParameter(6, 0.0 );
	
	
	
	BCLog::SetLogLevel(BCLog::debug);
	
	//BCHistogramFitter *histofitter = new BCHistogramFitter( ss_histoname.str().c_str(), tmphisto, ff_MCA );
	Gator::BcHistoFitterFast *histofitter = new Gator::BcHistoFitterFast( tmphisto, ff_MCA );
	
	// set options for MCMC
	//histofitter -> MCMCSetFlagPreRun (false);
	histofitter->MCMCSetPrecision(BCEngineMCMC::kHigh);
	//histofitter->SetFlagIntegration(true);
	//histofitter->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
	//histofitter->MCMCSetNIterationsPreRunMin(100000);
	histofitter->MCMCSetNIterationsRun(100000);
	//histofitter->MCMCSetMinimumEfficiency(0.05);
	//histofitter->MCMCSetMaximumEfficiency(0.2);
	
	
	int nChains = histofitter->MCMCGetNChains();//This depends on the precision
	
	//Used to generate random numbers between 0 and 1
	TRandom3 RdmGen(0);
	//set initial value of parameters
	vector<double> paramsval(nChains*nPars);
	for(int iChain = 0; iChain<nChains; iChain++){
		for(int iPar=0; iPar<nPars; iPar++){
			
			double val;
			BCParameter *Par = histofitter->GetParameter(iPar);
			
			double parmin, parmax;
			ff_MCA->GetParLimits(iPar, parmin, parmax);
			if(fixedpars.at(iPar) == true){
				val = ff_MCA->GetParameter(iPar);
				Par->SetLimits(val,val);
				Par->Fix(val);
				cout << "\nParameter <" << Par->GetName() << "> (fixed):" << endl;
				cout << "  Min = " << Par->GetLowerLimit() << endl;
				cout << "  Val = " << val << endl;
				cout << "  Max = " << Par->GetUpperLimit() << endl;
			}else{
				if(Par->Fixed()) Par->Unfix();
				ff_MCA->GetParLimits(iPar, parmin, parmax);
				val = parmin + RdmGen.Rndm()*(parmax-parmin);//Value between 0 and the parameter upper limit
				
				cout << "\nParameter <" << Par->GetName() << ">:" << endl;
				cout << "  Min = " << Par->GetLowerLimit() << endl;
				cout << "  Val = " << val << endl;
				cout << "  Max = " << Par->GetUpperLimit() << endl;
			}
			
			
			
			
			//This is for the case where the upper limit is infinite. r=0.5 will correspond to the initialized value.
			//double r = RdmGen.Rndm();//Value between 0 and 1
			//double val = ff_MCA->GetParameter(iPar)*r/(1-r);
			
			paramsval.at(iChain*nPars+iPar) = val;
		}
	}
	
	//exit(-1);
	
	histofitter->MCMCSetInitialPositions(paramsval);
	
	
	
	cout <<"\nMean init value: " << line.mean << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Mean")->GetLowerLimit() << " , " << histofitter->GetParameter("Mean")->GetUpperLimit() << ")" << endl;
	cout <<"\nAmpl init value: " << line.ampl << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Ampl")->GetLowerLimit() << " , " << histofitter->GetParameter("Ampl")->GetUpperLimit() << ")" << endl;
	cout <<"\nTail init value: " << line.tail << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Tail")->GetLowerLimit() << " , " << histofitter->GetParameter("Tail")->GetUpperLimit() << ")" << endl;
	cout <<"\nSigma init value: " << line.sigma << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Sigma")->GetLowerLimit() << " , " << histofitter->GetParameter("Sigma")->GetUpperLimit() << ")" << endl;
	cout <<"\nBeta init value: " << line.beta << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Beta")->GetLowerLimit() << " , " << histofitter->GetParameter("Beta")->GetUpperLimit() << ")" << endl;
	cout <<"\nStep init value: " << line.step << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Step")->GetLowerLimit() << " , " << histofitter->GetParameter("Step")->GetUpperLimit() << ")" << endl;
	cout <<"\nConstant init value: " << line.cost << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Constant")->GetLowerLimit() << " , " << histofitter->GetParameter("Constant")->GetUpperLimit() << ")" << endl;
	
	
	//exit(-1);
	
	histofitter->SetFlagIntegration(false);
	
	histofitter->Fit();
	
	histofitter->CalculatePValueLikelihood(histofitter->GetBestFitParameters());
	
	line.mean = histofitter->GetBestFitParameter(0);
	line.mean_err = histofitter->GetBestFitParameterError(0);
	line.ampl = histofitter->GetBestFitParameter(1);
	line.ampl_err = histofitter->GetBestFitParameterError(1);
	line.tail = histofitter->GetBestFitParameter(2);
	line.tail_err = histofitter->GetBestFitParameterError(2);
	//line.ratio = histofitter->GetBestFitParameter(2);
	//line.ratio_err = histofitter->GetBestFitParameterError(2);
	line.sigma = histofitter->GetBestFitParameter(3);
	line.sigma_err = histofitter->GetBestFitParameterError(3);
	line.beta = histofitter->GetBestFitParameter(4);
	line.beta_err = histofitter->GetBestFitParameterError(4);
	line.step = histofitter->GetBestFitParameter(5);
	line.step_err = histofitter->GetBestFitParameterError(5);
	line.cost = histofitter->GetBestFitParameter(6);
	line.cost_err = histofitter->GetBestFitParameterError(6);
	
	line.p_value = histofitter->GetPValue();
	line.p_value_ndof = histofitter->GetPValueNDoF();
	
	double nDof = (double)(histofitter->GetNDataPoints() - histofitter->GetNParameters());
	line.chi2ndof = ROOT::Math::chisquared_quantile_c( line.p_value_ndof, nDof );
	line.chi2 = line.chi2ndof*nDof;
	
	line.histo=tmphisto;
	line.fit=ff_MCA;
	
	
	
	//histofitter->DrawFit();
	TCanvas *c1 = new TCanvas("c1");
	c1->SetLogy();
	
	tmphisto->Draw();
	ff_MCA->Draw("same");
	
	c1->Update(); c1->Modified();
	
	if(theApp) theApp->Run(kTRUE);
	
	if(c1) delete c1;
	
	
	if(histofitter) delete histofitter;
	
	
	while(true){
		cout << "\nDo you want to keep this fit result for " << line.massN << "-" << line.element << " line at " << line.litEn << " keV? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
	}
	
	
	
	if(ans==string("y") || ans==string("Y")){
		return true;
	}else{
		
		while(true){
			cout << "\nDo you want to repeat the fit?\n" << "[y/n] > ";
			getline(cin,ans);
			if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
		}
		
		if(ans==string("y") || ans==string("Y")){
			cout << "\nCurrent parameters:\n" << endl;
			
			cout << "Mean         : " << line.mean << " +- " << line.mean_err << endl;
			cout << "Amplitude    : " << line.ampl << " +- " << line.ampl_err << endl;
			cout << "Tail (fixed) : " << line.tail << " +- " << line.tail_err << endl;
			cout << "Sigma        : " << line.sigma << " +- " << line.sigma_err << endl;
			cout << "Beta (fixed) : " << line.beta << " +- " << line.beta_err << endl;
			cout << "Step         : " << line.step << " +- " << line.step_err << endl;
			cout << "Constant     : " << line.cost << " +- " << line.cost_err << endl;
			
			
			if(!tail_flag){
				while(true){
					cout << "\nDo you want to reuse the current parameters?\n" << "[y/n] > ";
					getline(cin,ans);
					if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
				}
				
				if(ans==string("y") || ans==string("Y")){
					reuse_flag = true;
				}else{
					reuse_flag = false;
				}
			}
			
			//The parameters are kept only when the preliminary fit was performed without the tail part
			line.tail = intialValues.tail;
			line.beta = intialValues.beta;
			
			goto START;
		}
		
		return false;
	}
	
}


vector<CalibLine> LoadLinesFromTree(string datadir)
{
	
	TFile *treeFile = new TFile((datadir+string("fittedlines.root")).c_str(),"read");
	
	TTree *t1 = (TTree*)treeFile -> Get("linestree");
	
	vector<CalibLine> theLines;
	
	if(!t1){
		return theLines; //return empty vector
	}
	
	CalibLine tmp_line;
	
	//Variables to put in tree
	string *pmassN = &tmp_line.massN;
	string *pelement = &tmp_line.element;
	
	
	t1 -> SetBranchAddress("massNum",&pmassN);
	t1 -> SetBranchAddress("element",&pelement);
	t1 -> SetBranchAddress("litEn",&tmp_line.litEn);
	t1 -> SetBranchAddress("litEn_err",&tmp_line.litEn_err);
	t1 -> SetBranchAddress("mean",&tmp_line.mean);
	t1 -> SetBranchAddress("mean_err",&tmp_line.mean_err);
	t1 -> SetBranchAddress("sigma",&tmp_line.sigma);
	t1 -> SetBranchAddress("sigma_err",&tmp_line.sigma_err);
	t1 -> SetBranchAddress("beta",&tmp_line.beta);
	t1 -> SetBranchAddress("beta_err",&tmp_line.beta_err);
	t1 -> SetBranchAddress("ampl",&tmp_line.ampl);
	t1 -> SetBranchAddress("ampl_err",&tmp_line.ampl_err);
	t1 -> SetBranchAddress("tail",&tmp_line.tail);
	t1 -> SetBranchAddress("tail_err",&tmp_line.tail_err);
	//t1 -> SetBranchAddress("ratio",&tmp_line.ratio);
	//t1 -> SetBranchAddress("ratio_err",&tmp_line.ratio_err);
	t1 -> SetBranchAddress("step",&tmp_line.step);
	t1 -> SetBranchAddress("step_err",&tmp_line.step_err);
	t1 -> SetBranchAddress("cost",&tmp_line.cost);
	t1 -> SetBranchAddress("cost_err",&tmp_line.cost_err);
	t1 -> SetBranchAddress("p_value",&tmp_line.p_value);
	t1 -> SetBranchAddress("p_value_ndof",&tmp_line.p_value_ndof);
	
	
	int nLines = t1->GetEntries();
	
	for(int iLine=0; iLine<nLines; iLine++){
		t1 -> GetEntry(iLine);
		theLines.push_back(tmp_line);
	}
	
	
	return theLines;
}



TGraphErrors* plotResiduals(TF1* fit_pol2, TGraphErrors* grPlot, TFitResultPtr& resFit)
{
	
	int nPoints = grPlot->GetN();
	
	double* x = grPlot->GetX();
	double* ex = grPlot->GetEX();
	
	double* yTh = grPlot->GetY();
	double* eyTh = grPlot->GetEY();
	
	TGraphErrors* residual_gr = new TGraphErrors(nPoints);
	
	Double_t y, y_err, calc_val, calc_val_err;
	
	Double_t a, a_err, b, b_err, c, c_err, cov_ab, cov_ac, cov_bc;
	
	a = fit_pol2->GetParameter(0);
	a_err = fit_pol2->GetParError(0);
	b = fit_pol2->GetParameter(1);
	b_err = fit_pol2->GetParError(1);
	c = fit_pol2->GetParameter(2);
	c_err = fit_pol2->GetParError(2);
	
	cov_ab = resFit -> CovMatrix(0,1);
	cov_ac = resFit -> CovMatrix(0,2);
	cov_bc = resFit -> CovMatrix(1,2);
	
	
	for (Int_t i=0; i<nPoints; i++){
		
		calc_val = fit_pol2 -> Eval(x[i]);
		
		calc_val_err = TMath::Sqrt(a_err*a_err + b_err*x[i]*b_err*x[i] + ((c_err*x[i]*x[i]))*((c_err*x[i]*x[i])) + ((b+2*c*x[i])*ex[i])*((b+2*c*x[i])*ex[i]) + 2*(x[i])*(cov_ab + cov_ac*x[i] + cov_bc*x[i]*x[i]));
		
		y = calc_val - yTh[i];
		
		y_err = TMath::Sqrt((calc_val_err*calc_val_err)+(eyTh[i]*eyTh[i]));
		
		residual_gr -> SetPoint(i,x[i],y);
		residual_gr -> SetPointError(i,ex[i],y_err);
		
	}
	
	
	return residual_gr;
	
}


using Analysis::GatorCalib::GammaLineLikelihood;
using Analysis::Param;

bool doFitLL(TH1D* MCAhisto, CalibLine& line)
{
	GammaLineLikelihood *linelkh= new GammaLineLikelihood();
	
	//Make a copy of the ROI of the original histo and put it in a new histogram
	int firstbin = MCAhisto->FindBin(line.MCAlowch);
	int lastbin = MCAhisto->FindBin(line.MCAupch);
	int nbins = 1 + lastbin - firstbin;
	double xmin = MCAhisto->GetBinLowEdge(firstbin);
	double xmax = MCAhisto->GetBinLowEdge(lastbin+1);
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << line.element << line.massN << "_" << (int)(line.litEn+0.5) << "keV" ;
	TH1D* linehisto = new TH1D(ss_histoname.str().c_str(),";Channel;Counts",nbins,xmin,xmax);
	//linehisto -> SetDirectory(0);
	
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = MCAhisto->GetBinCenter(bin);
		linehisto->Fill(bincenter,MCAhisto->GetBinContent(bin));
	}
	
	line.histo = linehisto;
	
	TF1* ff_MCA = new TF1("ff_MCA", peakFitFuncA, xmin, xmax, 7);
	ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
	//ff_MCA->SetParNames("Mean","Ampl","Ratio","Sigma","Beta","Constant");
	
	line.fit = ff_MCA;
	
	//Initialization of critical parameters before the fitting
	amplInit(linehisto, line);
	costInit(linehisto, line);
	stepInit(linehisto, line);
	//line.step = line.ampl/1000;
	//line.step = 0.0;
	//sigmaInit(tmphisto,line);
	
	
	linelkh->Init(linehisto, &line);
	
	//linelkh->SetTrackVals();
	linelkh->SetEngineVerb(1);
	linelkh->SetVerbosity(2);
	//linelkh->SetMaxSteps(1000);
	linelkh->MetropolisMLE();
	linelkh->ScaleLogProbToMax();
	
	//This function copies also the results in the CalibLine structure
	linelkh->Minuit2MLE();
	
	
	ff_MCA -> SetParameter(0,line.mean);
	ff_MCA -> SetParameter(1,line.ampl); //Gaussian strenght
	//ff_MCA -> SetParameter(2,line.ampl*line.tail); //Tail strenght
	ff_MCA -> SetParameter(2,line.tail); //Tail ratio
	ff_MCA -> SetParameter(3,line.sigma);
	ff_MCA -> SetParameter(4,line.beta);
	ff_MCA -> SetParameter(5,line.step);
	ff_MCA -> SetParameter(6,line.cost);
	
	
	
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetLogy();
	
	linehisto->Draw();
	ff_MCA->Draw("SAME");
	
	c1->Update();
	
	if(theApp) theApp->Run(kTRUE);
	
	
	string ans("");
	while(true){
		cout << "\nDo you want to keep this fit result for " << line.massN << "-" << line.element << " line at " << line.litEn << " keV? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
	}
	
	
	
	bool retval;
	
	if(ans==string("y") || ans==string("Y")){
		retval = true;
	}else{
		retval = false;
	}
	
	if(linelkh){
		delete linelkh;
	}
	
	if(c1) delete c1;
	
	return retval;
}
