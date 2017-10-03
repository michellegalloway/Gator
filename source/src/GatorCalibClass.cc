#include "GatorCalibClass.hh"
//#include "loadSPE.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"

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

Gator::GatorCalib::GatorCalib()
{
	fHisto=NULL;
}


void Gator::GatorCalib::AddLine(CalibLine *line)
{
	if(line->linename == string(""))
	{
		stringstream ss_tmp;
		
		bool found_flag = false;
		unsigned counter=0;
		do{
			ss_tmp.str(""); ss_tmp << "Line" << fLinesmap.size()+counter; counter++;
		}while( fLinesmap.find(ss_tmp.str()) != fLinesmap.end() );
		
		fLinesmap[ss_tmp.str()] = line;
		return;
	}
	
	if(fLinesmap.find(line->linename) != fLinesmap.end()){
		if(fLinesmap[line->linename]) delete fLinesmap[line->linename];
	}
	fLinesmap[line->linename] = line;
}


void Gator::GatorCalib::LoadCalibFiles(const string& sourcename, const string& dir)
{
	if( fSpectramap.find(sourcename)!=fSpectramap.end())
	{
		if(fSpectramap[sourcename]) delete fSpectramap[sourcename];
	}
	
	double calibtime = 0.;
	cout << "Calibration dir: <" << dir << ">" << endl;
	TH1D* MCAhisto = loadSpe(dir.c_str(), calibtime);
	if(!MCAhisto) return;
	MCAhisto->SetTitle(";MCA channel; Counts");
	MCAhisto->SetStats(kFALSE);
	
	fSpectramap[sourcename] = MCAhisto;
}


TH1D* Gator::GatorCalib::SelectSpectrum(const string& sourcename)
{
	if( fSpectramap.find(sourcename)==fSpectramap.end()) fHisto = NULL; return fHisto;
	
	fHisto = fSpectramap[sourcename];
	
	return fHisto;
}


TH1D* Gator::GatorCalib::DrawSpectrum(const string& opt)
{
	if(fHisto) fHisto->Draw(opt.c_str());
	
	return fHisto;
}


TH1D* Gator::GatorCalib::GetSpectrum(const string& sourcename)
{
	if( fSpectramap.find(sourcename)==fSpectramap.end()) return NULL;
	
	return fSpectramap[sourcename];
}


Gator::CalibLine* Gator::GatorCalib::GetCalibLine(const string& linename)
{
	if( fLinesmap.find(linename)==fLinesmap.end()) return NULL;
	
	return fLinesmap[linename];
}


BCHistogramFitter* Gator::GatorCalib::FitLine(const string& linename)
{
	if( fLinesmap.find(linename)==fLinesmap.end()) return NULL;
	
	if(!fHisto) return NULL;
	
	CalibLine *line = fLinesmap[linename];
	
	
	//Make a copy of the ROI of the original histo and put it in a new histogram
	int firstbin = fHisto->FindBin(line->MCAlowch);
	int lastbin = fHisto->FindBin(line->MCAupch);
	int nbins = 1 + lastbin - firstbin;
	double xmin = fHisto->GetBinLowEdge(firstbin);
	double xmax = fHisto->GetBinLowEdge(lastbin+1);
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << (int)(line->litEn+0.5) << "keV" ;
	TH1D* tmphisto = new TH1D(ss_histoname.str().c_str(),";Channel;Counts",nbins,xmin,xmax);
	tmphisto -> SetDirectory(0);
	
	for(int bin=firstbin; bin <=lastbin; bin++){
		double bincenter = fHisto->GetBinCenter(bin);
		tmphisto->Fill(bincenter, fHisto->GetBinContent(bin));
	}
	
	
	//Start with the fitting process
	TF1* ff_MCA = new TF1("ff_MCA", &Gator::GatorCalib::peakFitFunc, xmin, xmax, 7);
	ff_MCA->SetParNames("Mean","Ampl","Tail","Sigma","Beta","Step","Constant");
	//ff_MCA->SetParNames("Mean","Ampl","Ratio","Sigma","Beta","Constant");
	
	costInit(tmphisto,*line);
	stepInit(tmphisto,*line);
	amplInit(tmphisto,*line);
	//sigmaInit(tmphisto,line);
	
	
	ff_MCA -> SetParLimits(0,line->MCAlowch,line->MCAupch);
	ff_MCA -> SetParLimits(1,0.,2*line->ampl); //Gaussian strenght
	//ff_MCA -> SetParLimits(2,0.,3*line->ampl); //Tail strenght
	//ff_MCA -> SetParLimits(2,0.,1.); //Tail ratio
	ff_MCA -> SetParLimits(3,0,2.5*line->sigma);
	ff_MCA -> SetParLimits(4,0.,10*line->beta);
	ff_MCA -> SetParLimits(5,0.,5*line->step);
	ff_MCA -> SetParLimits(6,0.,2*line->cost);
	
	
	
	BCLog::SetLogLevel(BCLog::debug);
	
	BCHistogramFitter *histofitter = new BCHistogramFitter(tmphisto, ff_MCA);
	
	// set options for MCMC
	//fHistofitter -> MCMCSetFlagPreRun (false);
	//fHistofitter->MCMCSetPrecision(BCEngineMCMC::kLow);
	//fHistofitter->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
	//fHistofitter->MCMCSetNIterationsPreRunMin(100000);
	histofitter->MCMCSetNIterationsRun(100000);
	
	
	//set initial value of parameters
	vector<double> paramsval;
	
	cout <<"\nMean init value: " << line->mean << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Mean")->GetLowerLimit() << " , " << histofitter->GetParameter("Mean")->GetUpperLimit() << ")" << endl;
	cout <<"\nAmpl init value: " << line->ampl << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Ampl")->GetLowerLimit() << " , " << histofitter->GetParameter("Ampl")->GetUpperLimit() << ")" << endl;
	cout <<"\nTail init value: " << line->tail << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Tail")->GetLowerLimit() << " , " << histofitter->GetParameter("Tail")->GetUpperLimit() << ")" << endl;
	cout <<"\nSigma init value: " << line->sigma << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Sigma")->GetLowerLimit() << " , " << histofitter->GetParameter("Sigma")->GetUpperLimit() << ")" << endl;
	cout <<"\nBeta init value: " << line->beta << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Beta")->GetLowerLimit() << " , " << histofitter->GetParameter("Beta")->GetUpperLimit() << ")" << endl;
	cout <<"\nStep init value: " << line->step << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Step")->GetLowerLimit() << " , " << histofitter->GetParameter("Step")->GetUpperLimit() << ")" << endl;
	cout <<"\nConstant init value: " << line->cost << endl;
	cout <<"\tLimits: (" << histofitter->GetParameter("Constant")->GetLowerLimit() << " , " << histofitter->GetParameter("Constant")->GetUpperLimit() << ")" << endl;
	
	int nPars = paramsval.size();
	int nChains = histofitter->MCMCGetNChains();//This depends on the precision
	
	for(int chain = 0; chain<nChains; chain++){
		paramsval.push_back(line->mean);
		paramsval.push_back(line->ampl);
		paramsval.push_back(line->tail);
		paramsval.push_back(line->sigma);
		paramsval.push_back(line->beta);
		paramsval.push_back(line->step);
		paramsval.push_back(line->cost);
	}
	histofitter->MCMCSetInitialPositions(paramsval);
	
	
	histofitter->Fit();
	
	line->mean = histofitter->GetBestFitParameter(0);
	line->mean_err = histofitter->GetBestFitParameterError(0);
	line->ampl = histofitter->GetBestFitParameter(1);
	line->ampl_err = histofitter->GetBestFitParameterError(1);
	line->tail = histofitter->GetBestFitParameter(2);
	line->tail_err = histofitter->GetBestFitParameterError(2);
	//line->ratio = histofitter->GetBestFitParameter(2);
	//line->ratio_err = histofitter->GetBestFitParameterError(2);
	line->sigma = histofitter->GetBestFitParameter(3);
	line->sigma_err = histofitter->GetBestFitParameterError(3);
	line->beta = histofitter->GetBestFitParameter(4);
	line->beta_err = histofitter->GetBestFitParameterError(4);
	line->step = histofitter->GetBestFitParameter(5);
	line->step_err = histofitter->GetBestFitParameterError(5);
	line->cost = histofitter->GetBestFitParameter(6);
	line->cost_err = histofitter->GetBestFitParameterError(6);
	line->p_value = histofitter->GetPValue();
	line->p_value_ndof = histofitter->GetPValueNDoF();
	
	histofitter->DrawFit("", true); // draw with a legend
	
	//Open a canvas to inspect the correlation of the beta and of the tail relative amplitude
	TPad *old_pad = (TPad*)gPad->cd();
	
	TCanvas *c2 = new TCanvas("c2");
	c2 -> cd();
	
	(histofitter->GetMarginalized("Tail","Beta"))->Draw();
	
	
	if(ff_MCA) delete ff_MCA;
	if(tmphisto) line->histo=tmphisto;
	if(c2) delete c2;
	
	old_pad->cd();
	
	return histofitter;
}


bool Gator::GatorCalib::LoadLinesFromTree(const string& rootfile)
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
		return false;
	}
	
	
	TTree* t1 = (TTree*)infile->Get("linestree");
	
	if(!t1)
	{
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


bool Gator::GatorCalib::SaveLines(const string& rootfile, bool update)
{
	if(access(rootfile.c_str(), W_OK) != 0)
	{
		return false;
	}
	
	if(!fLinesmap.size()) return false;
	
	
	TFile *outfile=NULL;
	
	if(update)
	{
		outfile->TFile::Open(rootfile.c_str(),"update");
	}else{
		outfile->TFile::Open(rootfile.c_str(),"recreate");
	}
	
	if( (!outfile) || (!outfile->IsWritable()) ) return false;
	
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
	
	t1 -> Branch("histo", "TH1D", &tmp_line.histo);
	
	
	map<string, CalibLine*>::iterator It;
	for(It=fLinesmap.begin(); It!=fLinesmap.end(); It++)
	{
		tmp_line = (*It->second);
		t1->Fill();
	}
	
	outfile->WriteTObject(t1);
	//delete t1;
	delete outfile;
	return true;
}


TH1D* Gator::GatorCalib::loadSpe(const char* dir, double& aqtime)
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






double Gator::GatorCalib::peakFitFunc(double* x, double* par){
	
	Double_t E,P,T,A,sigma,beta,S,C;
	
	E = x[0];
	
	P = par[0];
	A = par[1];
	T = par[2];
	sigma = par[3];
	beta = par[4];
	S = par[5];
	C = par[6];
	
	return A*( TMath::Exp( -pow((E-P)/sigma,2)/2 )/sigma + T*TMath::Exp( pow(sigma/beta,2)/2 + beta*(E-P) )*TMath::Erfc( ((E-P)*beta+pow(sigma,2) )/(sqrt(2)*sigma*beta) ) + S*TMath::Erfc((E-P)/(sqrt(2)*sigma)) ) + C;
	
	//return A*( TMath::Gaus(E,P,sigma) + R * TMath::Exp(beta*(E-P))*TMath::Erfc((E-(P-beta*sigma*sigma))/(sqrt(2)*sigma)) ) + C;
	
}
//////////////////////////////////////////////////////////////////////////////////


//----------Initialization functions for parameters-----------//

void Gator::GatorCalib::costInit(TH1D* histo, CalibLine& line){
	
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


void Gator::GatorCalib::stepInit(TH1D* histo, CalibLine& line){
	
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


void Gator::GatorCalib::amplInit(TH1D* histo, CalibLine& line){
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


void Gator::GatorCalib::sigmaInit(TH1D* histo, CalibLine& line){
	
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


#if defined(__CLING__)
//#include "loadSPE.cc"
#endif


