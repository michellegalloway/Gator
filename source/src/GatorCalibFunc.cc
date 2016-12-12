#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

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

#include "GatorCalibFitters.h"
#include "screenfncs.h"
#include "misc.h"

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

bool doFit(TH1D* MCAhisto, CalibLine& line)
{
	
	//Make a copy of the ROI of the original histo and put it in a new histogram
	int firstbin = MCAhisto->FindBin(line.MCAlowch);
	int lastbin = MCAhisto->FindBin(line.MCAupch);
	int nbins = 1 + lastbin - firstbin;
	double xmin = MCAhisto->GetBinLowEdge(firstbin);
	double xmax = MCAhisto->GetBinLowEdge(lastbin+1);
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << (int)(line.litEn+0.5) << "keV" ;
	TH1D* tmphisto = new TH1D(ss_histoname.str().c_str(),";Channel;Counts",nbins,xmin,xmax);
	tmphisto -> SetDirectory(0);
	
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = MCAhisto->GetBinCenter(bin);
		tmphisto->Fill(bincenter,MCAhisto->GetBinContent(bin));
	}
	
	TCanvas *c1 = new TCanvas("c1");
	c1->SetLogy();
	
	//Start with the fitting process
	TF1* ff_MCA = new TF1("ff_MCA",peakFitFunc,xmin,xmax,7);
	ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
	//ff_MCA->SetParNames("Mean","Ampl","Ratio","Sigma","Beta","Constant");
	
	costInit(tmphisto,line);
	stepInit(tmphisto,line);
	amplInit(tmphisto,line);
	//sigmaInit(tmphisto,line);
	
	
	ff_MCA -> SetParLimits(0,line.MCAlowch,line.MCAupch);
	ff_MCA -> SetParLimits(1,0.,2*line.ampl); //Gaussian strenght
	ff_MCA -> SetParLimits(2,0.,line.ampl/10); //Tail strenght
	//ff_MCA -> SetParLimits(2,0.,10.); //Tail ratio
	ff_MCA -> SetParLimits(3,0,2.5*line.sigma);
	ff_MCA -> SetParLimits(4,0.,10*line.beta);
	ff_MCA -> SetParLimits(5,0.,5*line.step);
	ff_MCA -> SetParLimits(6,0.,2*line.cost);
	
	
	
	BCLog::SetLogLevel(BCLog::debug);
	
	BCHistogramFitter *histofitter = new BCHistogramFitter(*((TH1*)tmphisto), *ff_MCA);
	
	// set options for MCMC
	//histofitter -> MCMCSetFlagPreRun (false);
	//histofitter->MCMCSetPrecision(BCEngineMCMC::kLow);
	//histofitter->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
	//histofitter->MCMCSetNIterationsPreRunMin(100000);
	//histofitter->MCMCSetNIterationsRun(100000);
	
	
	//set initial value of parameters
	vector<double> paramsval;
	
	cout <<"\nMean init value: " << line.mean << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Mean").GetLowerLimit() << " , " << histofitter->GetParameter("Mean").GetUpperLimit() << ")" << endl;
	
	cout <<"\nAmpl init value: " << line.ampl << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Ampl").GetLowerLimit() << " , " << histofitter->GetParameter("Ampl").GetUpperLimit() << ")" << endl;
	
	cout <<"\nTail init value: " << line.tail << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Tail").GetLowerLimit() << " , " << histofitter->GetParameter("Tail").GetUpperLimit() << ")" << endl;
	
	cout <<"\nSigma init value: " << line.sigma << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Sigma").GetLowerLimit() << " , " << histofitter->GetParameter("Sigma").GetUpperLimit() << ")" << endl;
	
	cout <<"\nBeta init value: " << line.beta << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Beta").GetLowerLimit() << " , " << histofitter->GetParameter("Beta").GetUpperLimit() << ")" << endl;
	
	cout <<"\nStep init value: " << line.step << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Step").GetLowerLimit() << " , " << histofitter->GetParameter("Step").GetUpperLimit() << ")" << endl;
	
	cout <<"\nConstant init value: " << line.cost << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Constant").GetLowerLimit() << " , " << histofitter->GetParameter("Constant").GetUpperLimit() << ")" << endl;
	
	
	int nChains = histofitter->GetNChains();//This depends on the precision
	
	for(int chain = 0; chain<nChains; chain++){
		paramsval.push_back(line.mean);
		paramsval.push_back(line.ampl);
		paramsval.push_back(line.tail);
		paramsval.push_back(line.sigma);
		paramsval.push_back(line.beta);
		paramsval.push_back(line.step);
		paramsval.push_back(line.cost);
	}
	int nPars = paramsval.size();
	
	
	histofitter->SetInitialPositions(paramsval);
	
	
	histofitter->Fit();
	
	line.mean = (histofitter->GetBestFitParameters()).at(0);
	line.mean_err = (histofitter->GetBestFitParameterErrors()).at(0);
	line.ampl = (histofitter->GetBestFitParameters()).at(1);
	line.ampl_err = (histofitter->GetBestFitParameterErrors()).at(1);
	line.tail = (histofitter->GetBestFitParameters()).at(2);
	line.tail_err = (histofitter->GetBestFitParameterErrors()).at(2);
	//line.ratio = histofitter->GetBestFitParameter(2);
	//line.ratio_err = histofitter->GetBestFitParameterError(2);
	line.sigma = (histofitter->GetBestFitParameters()).at(3);
	line.sigma_err = (histofitter->GetBestFitParameterErrors()).at(3);
	line.beta = (histofitter->GetBestFitParameters()).at(4);
	line.beta_err = (histofitter->GetBestFitParameterErrors()).at(4);
	line.step = (histofitter->GetBestFitParameters()).at(5);
	line.step_err = (histofitter->GetBestFitParameterErrors()).at(5);
	line.cost = (histofitter->GetBestFitParameters()).at(6);
	line.cost_err = (histofitter->GetBestFitParameterErrors()).at(6);
	line.p_value = histofitter->GetPValue();
	line.p_value_ndof = histofitter->GetPValueNDoF();
	
	histofitter->DrawFit("HIST", false); // draw with a legend
	
	c1->Update();
	
	if(theApp) theApp->Run(kTRUE);
	
	//Open a canvas to inspect the correlation of the beta and of the tail relative amplitude
	
	TCanvas *c2 = new TCanvas("c2");
	c2 -> cd();
	
	(histofitter->GetMarginalized("Tail","Beta")).Draw();
	
	//c2->Update();
	
	if(theApp) theApp->Run(kTRUE);
	
	
	string ans("");
	while(true){
		cout << "\nDo you want to keep this fit result for " << line.massN << "-" << line.element << " line at " << line.litEn << " keV? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
	}
	
	if(ff_MCA) delete ff_MCA;
	if(tmphisto) line.histo=tmphisto;
	if(c1) delete c1;
	if(c2) delete c2;
	
	if(ans==string("y") || ans==string("Y")){
		return true;
	}else{
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
