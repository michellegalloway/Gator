#include "GatorGlobals.hh"
#include "GatorCalibClass.hh"

#include "TCanvas.h"
#include "TApplication.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TPad.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"

#include "BAT/BCAux.h"
#include "BAT/BCLog.h"
#include "BAT/BCHistogramFitter.h"
#include "BAT/BCParameter.h"
#include "BAT/BCH2D.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>



GatorCalib::GatorCalib()
{
	fHisto=NULL;
	fDebug=false;
}

GatorCalib::~GatorCalib()
{
	return;
}

CalibLine* GatorCalib::AddLine(CalibLine *line, const string& spectrumname)
{
	if(fSpectramap.find(spectrumname) == fSpectramap.end()) return NULL;
	
	if(spectrumname == string("")) return NULL;
	
	TH1D* pSpect = fSpectramap[spectrumname];
	
	int firstbin = pSpect->FindBin( line->MCAlowch );
	int lastbin = pSpect->FindBin( line->MCAupch );
	int nbins = 1 + lastbin - firstbin;
	double xmin = pSpect->GetBinLowEdge(firstbin);
	double xmax = pSpect->GetBinLowEdge(lastbin+1);
	
	stringstream ss_tmp;
	ss_tmp.str(""); ss_tmp << line->massN << line->element << "_" << (int)(line->litEn+0.5) << "keV";
	
	line->linename = ss_tmp.str();
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << ((int)(line->litEn+0.5)) << "keV" ;
	line->histo = new TH1D(ss_histoname.str().c_str(),";Channel;Counts", nbins, xmin, xmax);
	line->histo->SetDirectory(0);
	
	for(int iBin=firstbin; iBin<=lastbin; iBin++)
	{
		line->histo->Fill( pSpect->GetBinCenter(iBin), pSpect->GetBinContent(iBin) );
	}
	
	line->fit = new TF1("ff_MCA", &GatorCalib::peakFitFuncB, xmin, xmax, 7);
	
	if(fLinesmap.find(line->linename) != fLinesmap.end()){
		if(fLinesmap[line->linename]) delete fLinesmap[line->linename];
	}
	fLinesmap[line->linename] = line;
	fLinesSpectraMap[line->linename] = spectrumname;
	
	return fLinesmap[line->linename];
}


CalibLine* GatorCalib::AddLine(const string& _massnum, const string& _element, const double& _litEn, const double& _litEnErr, const string& spectrum)
{
	if(fSpectramap.find(spectrum) == fSpectramap.end()) return NULL;
	
	if(spectrum == string("")) return NULL;
	
	stringstream linename; linename.str("");
	linename << _massnum << _element << "_" << ((int)(_litEn+0.5)) << "keV";
	return AddLine(new CalibLine(linename.str(), _massnum, _element, _litEn, _litEnErr), spectrum);
}


CalibLine* GatorCalib::AddLine(const string& _massnum, const string& _element, const double& _litEn, const string& spectrum)
{
	return AddLine(_massnum, _element, _litEn, 0.0, spectrum);
}


void GatorCalib::LoadCalibFiles(const string& sourcename, const string& dir)
{
	if( fSpectramap.find(sourcename)!=fSpectramap.end())
	{
		if(fSpectramap[sourcename]) delete fSpectramap[sourcename];
	}
	
	double calibtime = 0.;
	cout << "Calibration dir: <" << dir << ">" << endl;
	TH1D* MCAhisto = GatorCalib::loadSpe(dir.c_str(), calibtime);
	if(!MCAhisto) return;
	MCAhisto->SetName((sourcename+string("_MCAspec")).c_str());
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	fSpectramap[sourcename] = MCAhisto;
	fHisto = MCAhisto;
}


TH1D* GatorCalib::SelectSpectrum(const string& sourcename)
{
	if( fSpectramap.find(sourcename)==fSpectramap.end())
	{
		fHisto = NULL;
		return fHisto;
	}
	
	fHisto = fSpectramap[sourcename];
	
	return fHisto;
}


TH1D* GatorCalib::DrawSpectrum(const string& opt)
{
	if(fHisto) fHisto->Draw(opt.c_str());
	
	return fHisto;
}


TH1D* GatorCalib::GetSpectrum(const string& sourcename)
{
	if( fSpectramap.find(sourcename)==fSpectramap.end()) return NULL;
	
	return fSpectramap[sourcename];
}


CalibLine* GatorCalib::GetCalibLine(const string& linename)
{
	if( fLinesmap.find(linename)==fLinesmap.end()) return NULL;
	
	return fLinesmap[linename];
}


BCHistogramFitter* GatorCalib::FitLine(CalibLine& line)
{
	if(IsDebug())
	{
		cout << "\nDebug ---> Gator::GatorCalib::FitLine(...): entering the routine." << endl;
	}
	
	
	TH1D* tmphisto = line.histo;
	if(!tmphisto) return NULL;
	
	double xmin = tmphisto->GetBinLowEdge(1);
	double xmax = tmphisto->GetBinLowEdge(tmphisto->GetNbinsX()+1);
	
	
	TF1* ff_MCA = line.fit;
	if(!ff_MCA) return NULL;
	
	ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
	
	/*
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = MCAhisto->GetBinCenter(bin);
		tmphisto->Fill(bincenter,MCAhisto->GetBinContent(bin));
	}
	*/
	
	CalibLine intialValues = line;//Save the initial parameters as read from the file
	
	
	int nPars = ff_MCA->GetNpar();
	
	vector<bool> fixedpars(nPars,false);
	
	bool tail_flag = true;
	bool reuse_flag = false;
	
	START:
	
	
	if(IsDebug())
	{
		cout << "\nDebug ---> Gator::GatorCalib::FitLine(...): starting the fitting loop." << endl;
	}
	
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
	
	
	
	if(IsDebug())
	{
		cout << "\nDebug ---> Gator::GatorCalib::FitLine(...): setting the parameters init values in the fit function." << endl;
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
		if(IsDebug())
		{
			cout << "\nDebug ---> Gator::GatorCalib::FitLine(...): setting the parameters limits for the fit." << endl;
		}
		
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
	
	if(IsDebug())
	{
		cout << "\nDebug ---> Gator::GatorCalib::FitLine(...): instanciating the BCHistogramFitter object." << endl;
	}
	
	BcHistoFitterFast *histofitter = new BcHistoFitterFast( tmphisto, ff_MCA );
	
	// set options for MCMC
	//histofitter -> MCMCSetFlagPreRun (false);
	histofitter->MCMCSetPrecision(BCEngineMCMC::kVeryHigh);
	//histofitter->SetFlagIntegration(true);
	//histofitter->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
	//histofitter->MCMCSetNIterationsPreRunMin(100000);
	//histofitter->MCMCSetNIterationsRun(100000);
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
	tmphisto->Draw();
	ff_MCA->Draw("same");
	gPad->Update();
	
	//if(theApp) theApp->Run(kTRUE);
	
	
	while(true){
		cout << "\nDo you want to keep this fit result for " << line.massN << "-" << line.element << " line at " << line.litEn << " keV? [y/n]\n" << "> ";
		getline(cin,ans);
		if(ans!=string("y")  || ans!=string("Y") || ans!=string("n") || ans!=string("N")) break;
	}
	
	
	
	if(ans==string("y") || ans==string("Y")){
		return histofitter;
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
		
		return NULL;
	}
	
}



bool GatorCalib::LoadLinesFromTree(const string& rootfile)
{
	if(access(rootfile.c_str(), R_OK) != 0)
	{
		return false;
	}
	
	TFile *infile = TFile::Open(rootfile.c_str(), "read");
	
	if(!infile)
	{
		return false;
	}
	
	if(!infile->IsOpen())
	{
		delete infile;
		return false;
	}
	
	
	TTree* t1 = (TTree*)infile->Get("linestree");
	
	if(!t1)
	{
		delete infile;
		return false;
	}
	
	
	CalibLine tmp_line;
	
	t1 -> SetBranchAddress("massNum",&tmp_line.massN);
	t1 -> SetBranchAddress("element",&tmp_line.element);
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
	
	t1 -> SetBranchAddress("fit", &tmp_line.fit);
	t1 -> SetBranchAddress("histo", &tmp_line.histo);
	
	
	int nEvs = t1->GetEntries();
	
	stringstream ss_tmp;
	
	for(int iEv=0; iEv<nEvs; iEv++)
	{
		t1->GetEntry();
		
		ss_tmp.str("");
		ss_tmp << tmp_line.massN << tmp_line.element << (int)(tmp_line.litEn+0.5) << "keV";
		
		if( fLinesmap.find(ss_tmp.str())!=fLinesmap.end() )
		{
			if(fLinesmap[ss_tmp.str()]) delete fLinesmap[ss_tmp.str()];
		}
		
		fLinesmap[ss_tmp.str()] = new CalibLine(tmp_line);
	}
	
	t1->ResetBranchAddresses(); //This is to avoid eventual segmantation faults
	if(infile) delete infile; //This also closes the file and deletes the tree from the memory
	
	return true;
}


bool GatorCalib::SaveLines(const string& rootfile, bool update)
{
	if(access(rootfile.c_str(), W_OK) != 0)
	{
		cerr << "\nERROR---> Gator::GatorCalib::SaveLines(...): The file <" << rootfile << "> is not accessible!\n" << endl;
		return false;
	}
	
	if(!fLinesmap.size())
	{
		cerr << "\nERROR---> Gator::GatorCalib::SaveLines(...): The \"fLinesmap\" container is empty!\n" << endl;
		return false;
	}
	
	
	TFile *outfile=NULL;
	
	if(update)
	{
		outfile = TFile::Open(rootfile.c_str(),"update");
	}else{
		outfile = TFile::Open(rootfile.c_str(),"recreate");
	}
	
	if(!outfile)
	{
		cerr << "\nERROR---> Gator::GatorCalib::SaveLines(...): The \"outfile\" pointer in null!\n" << endl;
		return false;
	}
	
	if( !outfile->IsWritable() )
	{
		cerr << "\nERROR---> Gator::GatorCalib::SaveLines(...): The file <" << outfile->GetName() << "> is not writable!\n" << endl;
		return false;
	}
	
	TTree *t1 = new TTree("linestree","Results of the fits for each line");
	
	CalibLine tmp_line;
	
	t1 -> Branch("massNum","string",&tmp_line.massN);
	t1 -> Branch("element","string",&tmp_line.element);
	t1 -> Branch("litEn",&tmp_line.litEn,"litEn/D");
	t1 -> Branch("litEn_err",&tmp_line.litEn_err,"litEn_err/D");
	t1 -> Branch("mean",&tmp_line.mean,"mean/D");
	t1 -> Branch("mean_err",&tmp_line.mean_err,"mean_err/D");
	t1 -> Branch("sigma",&tmp_line.sigma,"sigma/D");
	t1 -> Branch("sigma_err",&tmp_line.sigma_err,"sigma_err/D");
	t1 -> Branch("beta",&tmp_line.beta,"beta/D");
	t1 -> Branch("beta_err",&tmp_line.beta_err,"beta_err/D");
	t1 -> Branch("ampl",&tmp_line.ampl,"ampl/D");
	t1 -> Branch("ampl_err",&tmp_line.ampl_err,"ampl_err/D");
	t1 -> Branch("tail",&tmp_line.tail,"tail/D");
	t1 -> Branch("tail_err",&tmp_line.tail_err,"tail_err/D");
	//t1 -> Branch("ratio",&tmp_line.ratio,"ratio/D");
	//t1 -> Branch("ratio_err",&tmp_line.ratio_err,"ratio_err/D");
	t1 -> Branch("step",&tmp_line.step,"step/D");
	t1 -> Branch("step_err",&tmp_line.step_err,"step_err/D");
	t1 -> Branch("cost",&tmp_line.cost,"cost/D");
	t1 -> Branch("cost_err",&tmp_line.cost_err,"cost_err/D");
	t1 -> Branch("p_value",&tmp_line.p_value,"p_value/D");
	t1 -> Branch("p_value_ndof",&tmp_line.p_value_ndof,"p_value_ndof/D");
	
	t1 -> Branch("histo", "TH1D", &tmp_line.histo, 32000, 0);
	t1 -> Branch("fit", "TF1", &tmp_line.fit, 32000, 0);
	
	map<string, CalibLine*>::iterator It;
	for(It=fLinesmap.begin(); It!=fLinesmap.end(); It++)
	{
		tmp_line = (*It->second);
		t1->Fill();
	}
	
	outfile->WriteTObject(t1, "linestree");
	//delete t1;
	delete outfile;
	return true;
}



TH1D* GatorCalib::loadSpe(const char* dir, double& aqtime)
{
	//-------------------------------------------------//
	// Load of the sample histogram from the SPE files //
	// this version is only for calibration files      //
	//-------------------------------------------------//


	string tmpstr("");
	stringstream tmpsstr;
	
	aqtime=0;
	
	Double_t time_tmp=0, entry=0; 
	Int_t channel=0;
	
	vector<Double_t>* entries_tmp = new vector<Double_t>;
	vector<Double_t>* entries = new vector<Double_t>;
	
	string dirname(dir);
	string command("ls ");
	vector<string> spelist;
	
	command += dirname;
	command += string("*.Spe > tmplist");
	
	system(command.c_str());
	ifstream tmplist("tmplist");
	
	while(getline(tmplist,tmpstr)){
		spelist.push_back(tmpstr);
	}
	
	int nSpe = spelist.size();
	
	system("rm -f tmplist");
	 
	ifstream spe;
	
	
	
	for(int iSpe=0; iSpe<nSpe; iSpe++){
		
		string infilename = spelist.at(iSpe);
		cout << "\nLoading file <" << infilename << ">" << endl;
		
		spe.open(infilename.c_str());
		
		if (spe.fail()){
			cout << "Error: can't open the file <" << infilename.c_str() << ">\nSkipping the entry!	" << endl;
			continue;//Skip the loop
		}
		
		entries_tmp->clear();
		
		int channels;
		
		for (Int_t line = 0; line<=11; line++) {
			//Here I go ahead in the file discarding the header lines
			getline(spe,tmpstr);
			
			if (line==9){//Acquisition live time (the first of the two numbers)
				//cout << "Content of \"tmpstr\" at line 9 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				time_tmp = atof(tmpstr.c_str());
				cout << "Acquisition time: " << time_tmp << " secs" << endl;
			}
			
			if (line==11){//Info on how many channels are acquired (the second of the two numbers)
				//cout << "Content of \"tmpstr\" at line 11 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				tmpsstr >> tmpstr;
				channels = atoi(tmpstr.c_str());
				channels++;//The number of channels is 1 more of the last channel number
				cout << "Number of channels: " << channels << endl;
			}
		}
		
		//Loading the ADC counts from the file and put them in the vector entry_tmp
		while (spe.good()) {
			spe >> tmpstr;
			entry = atof(tmpstr.c_str());
			if(iSpe==0){
				entries -> push_back(entry);
				if (entries->size()==channels) break;
			} else {
				entries_tmp -> push_back(entry);
				if (entries_tmp->size()==channels) break;
			}
				
			if (!spe.good()) break;
		}
		
		if(iSpe==0){
			cout << "Data vector size: " << entries->size() << endl;
		}else{
			cout << "Data vector size: " << entries_tmp->size() << endl;
			for(Int_t i=0; i<entries_tmp->size(); i++){
				(*entries)[i] += entries_tmp->at(i);
			}
		}
		
		spe.close();
		
		aqtime += time_tmp; //Stored total total acquisition time
		
	}
	
	//Create and fill the histogram
	TH1D* histoMCA = new TH1D("histoMCA","Spectrum ADC",entries->size(),1,entries_tmp->size());
	for(Int_t i=0; i<entries->size(); i++){
		histoMCA->Fill(i+1,entries->at(i));
	}
	
	//-------------------------------------------------------------//	
	// Finished to load of the sample histogram from the SPE files //
	//-------------------------------------------------------------//

	
	delete entries_tmp;
	delete entries;
	
	return histoMCA;
}






double GatorCalib::peakFitFuncA(double* x, double* par){
	
	double E,P,T,A,sigma,beta,S,C;
	
	E = x[0];
	
	P = par[0];
	A = par[1];
	T = par[2];
	sigma = par[3];
	beta = par[4];
	S = par[5];
	C = par[6];
	
	double gauss = TMath::Exp( -pow((E-P)/sigma,2)/2 );
	//double tail = TMath::Exp( (pow(sigma/beta,2)/2) + ((E-P)/beta) ) * TMath::Erfc(  ( (E-P)*beta+pow(sigma,2) )/( sqrt(2)*sigma*beta )  );
	//double tail = TMath::Exp( ((E-P)/beta) ) * TMath::Erfc(  ( (E-P)*beta+pow(sigma,2) )/( sqrt(2)*sigma*beta )  );
	double tail = TMath::Exp( ((E-P)/beta) ) * TMath::Erfc( (((E-P)/sigma) + (sigma/beta))/sqrt(2) );
	double step = TMath::Erfc((E-P)/(sqrt(2)*sigma));
	
	return A*(gauss + T*tail) + S*step + C;
	
	//return A*( TMath::Gaus(E,P,sigma) + R * TMath::Exp(beta*(E-P))*TMath::Erfc((E-(P-beta*sigma*sigma))/(sqrt(2)*sigma)) ) + C;
}
//////////////////////////////////////////////////////////////////////////////////


double GatorCalib::peakFitFuncB(double* x, double* par){
	
	double E,P,T,A,sigma,Gamma,S,C;
	
	E = x[0];
	
	P = par[0];
	A = par[1];
	T = par[2];
	sigma = par[3];
	Gamma = par[4];
	S = par[5];
	C = par[6];
	
	double gauss = TMath::Gaus( E,P,sigma );
	//double tail = TMath::Exp( (pow(sigma/beta,2)/2) + ((E-P)/beta) ) * TMath::Erfc(  ( (E-P)*beta+pow(sigma,2) )/( sqrt(2)*sigma*beta )  );
	//double tail = TMath::Exp( ((E-P)/beta) ) * TMath::Erfc(  ( (E-P)*beta+pow(sigma,2) )/( sqrt(2)*sigma*beta )  );
	double tail = TMath::Exp( ((E-P)*Gamma) ) * TMath::Erfc( ( E-(P-Gamma*sigma*sigma) )/(sqrt(2)*sigma) );
	double step = TMath::Erfc( (E-P)/(sqrt(2)*sigma) );
	
	//return A*((1-T)*gauss + T*tail) + S*step + C;
	return A*(gauss + T*tail) + S*step + C;
	
	//return A*( TMath::Gaus(E,P,sigma) + R * TMath::Exp(beta*(E-P))*TMath::Erfc((E-(P-beta*sigma*sigma))/(sqrt(2)*sigma)) ) + C;
	
}



//----------Initialization functions for parameters-----------//

void GatorCalib::costInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	
	Double_t left_mean = 0.;
	Double_t right_mean = 0.;
	
	for(int bin=firstbin; bin<firstbin+10; bin++){
		left_mean += (histo->GetBinContent(bin))/10.;
	}
	for(int bin=lastbin-9; bin<=lastbin; bin++){
		right_mean += (histo->GetBinContent(bin))/10.;
	}
	
	line.cost = TMath::Min(left_mean,right_mean);
	
	return;
}
//////////////////////////////////////////////////////////////////////////////////////


void GatorCalib::stepInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	
	Double_t left_mean = 0.;
	Double_t right_mean = 0.;
	
	for(int bin=firstbin; bin<firstbin+10; bin++){
		left_mean += (histo->GetBinContent(bin))/10.;
	}
	for(int bin=lastbin-9; bin<=lastbin; bin++){
		right_mean += (histo->GetBinContent(bin))/10.;
	}
	
	if((left_mean - right_mean)<=0){
		line.step = right_mean;
	}else{
		line.step = (left_mean - right_mean);
	}
	
	return;
}
/////////////////////////////////////////////////////////////////////////////////////


void GatorCalib::amplInit(TH1D* histo, CalibLine& line){
	//This should always be run after the function "costInit"
	double max =0.;
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	for(int bin=firstbin; bin<=lastbin; bin++){
		if(max < histo->GetBinContent(bin)) max = histo->GetBinContent(bin);
	}
	
	line.ampl = max-line.cost;
	//line.tail = 0.1*line.ampl;
	
	return;
}
//////////////////////////////////////////////////////////////////////////////////////


void GatorCalib::sigmaInit(TH1D* histo, CalibLine& line){
	
	Double_t firstbin = histo->FindBin(line.MCAlowch); 
	Double_t lastbin = histo->FindBin(line.MCAupch);
	
	double max =0.;
	int maxbin = firstbin;
	
	for(int bin=firstbin; bin<=lastbin; bin++){
		if(max < histo->GetBinContent(bin)){
			max = histo->GetBinContent(bin);
			maxbin = bin;
		}
	}
	
	int bin = maxbin;
	while(histo->GetBinContent(bin)>max/2.){
		bin--;
	}
	double lowval = histo->GetBinCenter(bin);
	
	bin = maxbin;
	while(histo->GetBinContent(bin)>max/2.){
		bin++;
	}
	double upval = histo->GetBinCenter(bin);
		
	line.sigma = (upval - lowval)/2; //This is not the sigma but a good starting point for the fit
	
	return;
}
/////////////////////////////////////////////////////////////////////////////////////



