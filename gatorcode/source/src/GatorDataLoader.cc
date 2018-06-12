#include "screenfncs.h"
#include "GatorStructs.h"
#include "GatorCounter.hh"

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TH1D.h"

#include "GatorDataLoader.hh"

using namespace std;

//Default Constructor
GatorDataLoader::GatorDataLoader(){
	
	md_aqtime = 0;
	md_Bgtime = 0;
	
	m_datadir=string("");
	m_datacalibdir=string("");
	
	mp_data_calib_bw = NULL;
	mp_data_calib_fw = NULL;
	mp_data_MCA_resol = NULL;
	
	m_Bgdir=string("");
	m_Bgcalibdir=string("");
	
	mp_Bgcalib_bw = NULL;
	mp_Bgcalib_fw = NULL;
	mp_BgMCA_resol = NULL;
	
	
	md_WidthScale = 3.0;
	
	mb_dataloaded = false;
	mb_Bgloaded = false;
	mb_linesloaded = false;
		
}


bool GatorDataLoader::LoadLines(string cfgfilename){
	
	//cout << "\nEntering in GatorDataLoader::LoadLine\n" << endl;
	
	if(!fexist(cfgfilename)){
		cerr << "GatorDataLoader::LoadLines --> ERROR::: The file <" << cfgfilename << "> doesn't exists!" << endl;
		mb_linesloaded = false;
		return mb_linesloaded;
	}
	
	if(mv_lines->size()>0) mv_lines->clear();
	
	LineStruct tmpLine;
	
	ifstream cfgFile(cfgfilename.c_str());
	string tmpstr;
	stringstream tmpss;
	while( getline(cfgFile,tmpstr) ){
		
		//cout << "\n\n" << tmpstr << endl;
		
		tmpss.clear();
		tmpss.str("");
		tmpss << tmpstr;
		
		tmpss >> tmpstr;//The Mass number
		tmpLine.mass = tmpstr;
		tmpss >> tmpstr;//The element name
		tmpLine.element = tmpstr;
		tmpss >> tmpstr;//The litterature energy
		tmpLine.Energy = strtod(tmpstr.c_str(),NULL);
		tmpss >> tmpstr;//The Branching Ratio
		tmpLine.BR = strtod(tmpstr.c_str(),NULL);
		tmpss >> tmpstr;//The Line Efficiency
		tmpLine.Effic = strtod(tmpstr.c_str(),NULL);
		tmpss >> tmpstr;//The product of the Line Efficiency and the BR (what is actually used)
		tmpLine.BRxEffic = strtod(tmpstr.c_str(),NULL);
		
		
		DefaultIntervals(tmpLine);
		
		mv_lines->push_back(tmpLine);
	}
	cfgFile.close();
	
	mb_linesloaded = mv_lines->size()>0;
	
	//cout << "\nExiting from GatorDataLoader::LoadLine\n" << endl;
	
	return mb_linesloaded;
}


void GatorDataLoader::DefaultIntervals(){
	
	int nLines = GetNLines();
	
	if(nLines==0) return;
	
	for(int iLine=0; iLine<nLines; iLine++){
		
		DefaultIntervals( *(GetLine(iLine)) );
		
	}
	
	return;
	
}


void GatorDataLoader::DefaultIntervals(LineStruct &line){
	
	//cout << "\nEntering in GatorDataLoader::DefaultIntervals(LineStruct &line)\n" << endl;
	
	if(!mb_dataloaded) return;
	
	line.MCAcenter = mp_data_calib_bw->Eval(line.Energy);
	line.MCAsigma = mp_data_MCA_resol->Eval(line.MCAcenter);
	
	line.LowEdgeL = line.MCAcenter - 2*md_WidthScale*line.MCAsigma;
	line.UpEdgeL = line.MCAcenter - md_WidthScale*line.MCAsigma;
	line.LowEdgeS = line.UpEdgeL;
	line.UpEdgeS = line.MCAcenter + md_WidthScale*line.MCAsigma;
	line.LowEdgeR = line.UpEdgeS;
	line.UpEdgeR = line.MCAcenter + 2*md_WidthScale*line.MCAsigma;
	
	if(mb_BGflag){
		if(!mb_Bgloaded) return;
		
		line.BgMCAcenter = mp_Bgcalib_bw->Eval(line.Energy);
		line.BgMCAsigma = mp_BgMCA_resol->Eval(line.BgMCAcenter);
	
		line.BgLowEdgeL = line.BgMCAcenter - 2*md_WidthScale*line.BgMCAsigma;
		line.BgUpEdgeL = line.BgMCAcenter - md_WidthScale*line.BgMCAsigma;
		line.BgLowEdgeS = line.BgUpEdgeL;
		line.BgUpEdgeS = line.BgMCAcenter + md_WidthScale*line.BgMCAsigma;
		line.BgLowEdgeR = line.BgUpEdgeS;
		line.BgUpEdgeR = line.BgMCAcenter + 2*md_WidthScale*line.BgMCAsigma;
		
	}
	
	//cout << "  " << line.MCAcenter << endl;
	
	//cout << "\nExiting from GatorDataLoader::DefaultIntervals(LineStruct &line)\n" << endl;
	
	return;
}


void GatorDataLoader::SetLeftCompRegion(double lowEn, double upEn, int iLine){
	
	if(!mb_linesloaded) return;
	
	if(!mb_dataloaded) return;
	
	int nLines = GetNLines();
	
	if( (nLines<=0)||(iLine>=nLines) ) return;
	
	LineStruct *line = GetLine(iLine);
		
	
	line->LowEdgeL = mp_data_calib_bw->Eval(lowEn);
	line->UpEdgeL = mp_data_calib_bw->Eval(upEn);
	
	if(mb_BGflag){
		if(!mb_Bgloaded) return;
		
		line->BgLowEdgeL = mp_Bgcalib_bw->Eval(lowEn);
		line->BgUpEdgeL = mp_Bgcalib_bw->Eval(upEn);
		
	}
	
}


void GatorDataLoader::SetSignalRegion(double lowEn, double upEn, int iLine){
	
	if(!mb_linesloaded) return;
	
	if(!mb_dataloaded) return;
	
	int nLines = GetNLines();
	
	if( (nLines<=0)||(iLine>=nLines) ) return;
	
	LineStruct *line = GetLine(iLine);
		
	
	line->LowEdgeS = mp_data_calib_bw->Eval(lowEn);
	line->UpEdgeS = mp_data_calib_bw->Eval(upEn);
	
	if(mb_BGflag){
		if(!mb_Bgloaded) return;
		
		line->BgLowEdgeS = mp_Bgcalib_bw->Eval(lowEn);
		line->BgUpEdgeS = mp_Bgcalib_bw->Eval(upEn);
		
	}
	
}


void GatorDataLoader::SetRightCompRegion(double lowEn, double upEn, int iLine){
	
	if(!mb_linesloaded) return;
	
	if(!mb_dataloaded) return;
	
	int nLines = GetNLines();
	
	if( (nLines<=0)||(iLine>=nLines) ) return;
	
	LineStruct *line = GetLine(iLine);
		
	
	line->LowEdgeR = mp_data_calib_bw->Eval(lowEn);
	line->UpEdgeR = mp_data_calib_bw->Eval(upEn);
	
	if(mb_BGflag){
		if(!mb_Bgloaded) return;
		
		line->BgLowEdgeR = mp_Bgcalib_bw->Eval(lowEn);
		line->BgUpEdgeR = mp_Bgcalib_bw->Eval(upEn);
		
	}
	
}


void GatorDataLoader::SetSampleDataDir(string datadir){
	
	if(datadir.substr(datadir.size()-1,1)!=string("/")) datadir = datadir + string("/");
	
	m_datadir = datadir;
	
	return;
}


bool GatorDataLoader::LoadData(string datadir){
	
	
	cout << datadir << endl;
	
	if(datadir!=string("")){
		if(datadir.substr(datadir.size()-1,1)!=string("/")) datadir = datadir + string("/");
		SetSampleDataDir(datadir);
	}
	if(m_datadir==string("")){
		cerr << "GatorDataLoader::LoadData --> ERROR::: No data directory provided!" << endl;
		mb_dataloaded=false;
		return mb_dataloaded;
	}
	
	//This part doesn't work very well!
	if(fexist(datadir+string(".calib"))){
		
		ifstream calibrepository( (datadir+string(".calibdir")).c_str() );
		getline(calibrepository, m_datacalibdir);
		if(m_datacalibdir.substr(m_datacalibdir.size()-1,1)!=string("/")) m_datacalibdir = m_datacalibdir + string("/");
		calibrepository.close();
		
	}else if(m_datacalibdir==string("")){
		
		cerr << "GatorDataLoader::LoadData --> ERROR::: No calibration directory provided!" << endl;
		mb_dataloaded=false;
		return mb_dataloaded;
		
	}
	
	mp_sampleSpect = loadSPE(m_datadir.c_str(), md_aqtime);
	mp_sampleSpect -> SetName("SampleMCAspec");
	mp_sampleSpect -> SetTitle("; MCA channel; counts");
	
	mb_dataloaded = SetSampleCalib(m_datacalibdir);
	
	return mb_dataloaded;
}


bool GatorDataLoader::SetSampleCalib(string calibdir){
	
	if(calibdir==string("")){
		cerr << "GatorDataLoader::SetSampleCalib --> ERROR::: The directory for the data calibration file cannot be an empty string!" << endl;
		return false;
	}
	
	if(calibdir.substr(calibdir.size()-1,1)!=string("/")) calibdir = calibdir + string("/");
	
	if(!(fexist(calibdir+string("calibration.root")) && fexist(calibdir+string("resolution.root")))){
		cerr << "GatorDataLoader::SetSampleCalib --> ERROR::: At least one of these files doesn't exist:\n  <" << calibdir+string("calibration.root") << ">\n  <" << calibdir+string("resolution.root") << ">" << endl;
		return false;
	}
	
	m_datacalibdir = calibdir;	
	
	TFile *calibfile = new TFile((m_datacalibdir+string("calibration.root")).c_str(), "read");
	mp_data_calib_bw = (TF1*)calibfile -> Get("calib_fcn_bw");
	mp_data_calib_fw = (TF1*)calibfile -> Get("calib_fcn_fw");
	calibfile->Close();
	TFile *resolfile = new TFile((m_datacalibdir+string("resolution.root")).c_str(), "read");
	mp_data_MCA_resol = (TF1*)resolfile -> Get("resol_ADC_func");
	resolfile->Close();
	
	if(calibfile) delete calibfile;
	if(resolfile) delete resolfile;
	
	return( mp_data_calib_bw && mp_data_calib_fw && mp_data_MCA_resol );
}


void GatorDataLoader::SetBackgroundDataDir(string Bgdir){
	
	if(Bgdir.substr(Bgdir.size()-1,1)!=string("/")) Bgdir = Bgdir + string("/");
	
	m_Bgdir = Bgdir;
	
	return;
}


bool GatorDataLoader::LoadBackground(string Bgdir){
	
	if(Bgdir!=string("")){
		if(Bgdir.substr(Bgdir.size()-1,1)!=string("/")) Bgdir = Bgdir + string("/");
		SetBackgroundDataDir(Bgdir);
	}
	if(m_Bgdir==string("")){
		cerr << "GatorDataLoader::LoadBackground --> ERROR::: No data directory provided!" << endl;
		mb_Bgloaded=false;
		return mb_Bgloaded;
	}
	
	if(fexist(Bgdir+string(".calib"))){
		
		ifstream calibrepository( (Bgdir+string(".calibdir")).c_str() );
		getline(calibrepository, m_Bgcalibdir);
		if(m_Bgcalibdir.substr(m_Bgcalibdir.size()-1,1)!=string("/")) m_Bgcalibdir = m_Bgcalibdir + string("/");
		calibrepository.close();
		
	}else if(m_Bgcalibdir ==string("")){
		
		cerr << "GatorDataLoader::SetBackgroundDataDir --> ERROR::: No calibration directory provided for background!" << endl;
		mb_Bgloaded=false;
		return mb_Bgloaded;
		
	}
	
	mp_backgroundSpect = loadSPE(Bgdir.c_str(), md_Bgtime);
	mp_backgroundSpect -> SetName("BgMCAspec");
	mp_backgroundSpect -> SetTitle("; MCA channel; counts");
	
	mb_Bgloaded = SetBackgroundCalib(m_Bgcalibdir);
	
	SetBackgroundFlag(mb_Bgloaded);
	
	return mb_Bgloaded;
}


bool GatorDataLoader::SetBackgroundCalib(string calibdir){
	
	if(calibdir==string("")){
		cerr << "GatorDataLoader::SetBackgroundCalib --> ERROR::: The directory for the background calibration file cannot be an empty string!" << endl;
		return false;
	}
	
	if(calibdir.substr(calibdir.size()-1,1)!=string("/")) calibdir = calibdir + string("/");
	
	if( !(fexist(calibdir+string("calibration.root")) && fexist(calibdir+string("resolution.root"))) ){
		cerr << "GatorDataLoader::SetBackgroundCalib --> ERROR::: At least one of these files doesn't exist:\n  <" << calibdir+string("calibration.root") << ">\n  <" << calibdir+string("resolution.root") << ">" << endl;
		return false;
	}
	
	m_Bgcalibdir = calibdir;
	
	TFile *calibfile = new TFile((m_Bgcalibdir+string("calibration.root")).c_str(), "read");
	mp_Bgcalib_bw = (TF1*)calibfile -> Get("calib_fcn_bw");
	mp_Bgcalib_fw = (TF1*)calibfile -> Get("calib_fcn_fw");
	calibfile->Close();
	TFile *resolfile = new TFile((m_Bgcalibdir+string("resolution.root")).c_str(), "read");
	mp_BgMCA_resol = (TF1*)resolfile -> Get("resol_ADC_func");
	resolfile->Close();
	
	if(calibfile) delete calibfile;
	if(resolfile) delete resolfile;
	
	return( mp_Bgcalib_bw && mp_Bgcalib_fw && mp_BgMCA_resol );
	
}





