#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"

#include "GatorSampleAnalysis.hh"
#include "screenfncs.h"


using namespace std;

//Default Constructor
GatorSampleAnalysis::GatorSampleAnalysis(string workdir){
	
	mb_dataAnalysed = false;
	
	m_workdir = workdir;
	
	m_logfilename = m_workdir+string("SampleAnalysis.log");
	
	md_quant = 0;
	
	m_unit = string("kg");
	
	mp_sampleEnSpect = NULL;
	mp_backgroundEnSpect = NULL;
	gr_DetActiv = NULL;
	gr_ULActiv = NULL;
}


void GatorSampleAnalysis::doAnalysis(){
	
	
	if(md_quant<0){
		cerr << "GatorSampleAnalysis::doAnalysis --> ERROR::: Mass value not allowed!" << endl;
		return;
	}
	
	if(!mb_linesloaded){
		cerr << "GatorSampleAnalysis::doAnalysis --> ERROR::: No lines loaded, impossible to perform any analysys!" << endl;
		return;
	}
	
	if( !(mb_dataloaded) ){
		cerr << "GatorSampleAnalysis::doAnalysis --> ERROR::: Sample data not loaded!" << endl;
		return;
	}
	
	if(mb_BGflag){
		if( !(mb_Bgloaded) ){
			cerr << "GatorSampleAnalysis::doAnalysis --> ERROR::: Background data not loaded!" << endl;
			return;
		}
	}
	
	int nLines = GetNLines();
	
	ofstream logfile(m_logfilename.c_str());
	
	logfile << "Number of lines to be analyzed: " << nLines << endl << endl;
	
	
	for(int iLine=0; iLine<nLines; iLine++){
		LineStruct *line = GetLine(iLine);
		
		LineCounts(*line);
		
		double BRxEff = line->BRxEffic;
		
		//Start from the background (check if there is a peak)
		if(mb_BGflag){
			
			double bg_line_interv = line->BgbinS;
			double bg_comp_interv = line->BgbinL+line->BgbinR;
			
			line->BgPeakCnts = line->BgcountsS;
			line->BgPeakCnts_err = sqrt(line->BgcountsS);
		
			line->BgCompCnts = (line->BgcountsL+line->BgcountsR)*bg_line_interv/bg_comp_interv;
			line->BgCompCnts_err = sqrt(line->BgcountsL+line->BgcountsR)*bg_line_interv/bg_comp_interv;
			
			line->BgLineCnts = line->BgPeakCnts - line->BgCompCnts;
			line->BgLineCnts_err = sqrt( pow(line->BgPeakCnts_err,2) + pow(line->BgCompCnts_err,2) );
			
			line->BgLdCnts = 2.86 + 4.87*TMath::Sqrt(1.36 + line->BgCompCnts); //From the Gator paper;
			
			/*
			cout << "\n\n\n" << endl;
			cout << bg_line_interv << endl;
			cout << bg_comp_interv << endl;
			cout << line->BgPeakCnts << endl;
			cout << line->BgCompCnts << endl;
			cout << line->BgLineCnts << endl;
			cout << "\n\n\n" << endl;
			*/
			
			if(line->BgLdCnts<line->BgLineCnts){
				line->BgDetected=true;
			}else{
				line->BgDetected=false;
			}
		}else{
			line->BgDetected=false;
		}
		
		
		double line_interv = line->binS;
		double comp_interv = line->binL+line->binR;
		
		
		line->PeakCnts=line->countsS;
		line->PeakCnts_err=sqrt(line->countsS);
		
		line->CompCnts = (line->countsL+line->countsR)*line_interv/comp_interv;
		line->CompCnts_err = sqrt(line->countsL+line->countsR)*line_interv/comp_interv;
		
		if(line->BgDetected){
			double bgcnts = (md_aqtime/md_Bgtime)*line->BgLineCnts;
			double bgcnts_err = (md_aqtime/md_Bgtime)*line->BgLineCnts_err;
			
			line->LineCnts = line->PeakCnts - line->CompCnts - bgcnts;
			line->LineCnts_err = sqrt( line->PeakCnts + line->CompCnts + pow(bgcnts_err,2) );
			
			line->LdCnts = 2.86 + 4.87*TMath::Sqrt(1.36 + line->CompCnts + bgcnts); //From the Gator paper
			
		}else{
			line->LineCnts = line->PeakCnts - line->CompCnts;
			line->LineCnts_err = sqrt( line->PeakCnts + line->CompCnts );
			
			line->LdCnts = 2.86 + 4.87*TMath::Sqrt(1.36 + line->CompCnts); //From the Gator paper
		}
		
		line->LdActiv = 1000*line->LdCnts/(md_quant*BRxEff*md_aqtime);
		
		line->Activ = 1000*line->LineCnts/(md_quant*BRxEff*md_aqtime);
		line->Activ_err = line->Activ*sqrt( 0.01 + pow(line->LineCnts_err/line->LineCnts,2) ); //Adding a 10% sistematic from simulations (From Gator Paper)
		
		if(line->LdCnts < line->LineCnts){
			line->Detected=true;
		}else{
			line->Detected=false;
		}
		
		
		
		logfile << line->mass << line->element << " line at " << line->Energy << " keV:" << endl;
		
		if(mb_BGflag){
			logfile << "Background:" << endl;
			logfile << " MCA center = " << line->BgMCAcenter << endl;
			logfile << " MCA sigma = " << line->BgMCAsigma << endl;
			logfile << " Left compton region = (" << line->BgLowEdgeL << "," << line->BgUpEdgeL << ") --> " << line->BgbinL << " chs" << endl;
			logfile << " Left compton counts = " << line->BgcountsL << endl;
			logfile << " Signal region = (" << line->BgLowEdgeS << "," << line->BgUpEdgeS << ") --> " << line->BgbinS << " chs" << endl;
			logfile << " Signal region counts = " << line->BgcountsS << endl;
			logfile << " Right compton region = (" << line->BgLowEdgeR << "," << line->BgUpEdgeR << ") --> " << line->BgbinR << " chs" << endl;
			logfile << " Right compton counts = " << line->BgcountsR << endl;
			logfile << " Estimated Compton counts in signal region = " << line->BgCompCnts << " +- " << line->BgCompCnts_err << endl;
			logfile << " Estimated line net counts = " << line->BgLineCnts << " +- " << line->BgLineCnts_err << endl;
			logfile << " Background detection limit = " << line->BgLdCnts << endl;
			if(line->BgDetected) logfile << " Background detected! " << endl;
			logfile << "------------" << endl ;
		}
		
		logfile << "Sample:" << endl;
		logfile << " MCA center = " << line->MCAcenter << endl;
		logfile << " MCA sigma = " << line->MCAsigma << endl;
		logfile << " Left compton region = (" << line->LowEdgeL << "," << line->UpEdgeL << ") --> " << line->binL << " chs" << endl;
		logfile << " Left compton counts = " << line->countsL << endl;
		logfile << " Signal region = (" << line->LowEdgeS << "," << line->UpEdgeS << ") --> " << line->binS << " chs" << endl;
		logfile << " Signal region counts = " << line->countsS << endl;
		logfile << " Right compton region = (" << line->LowEdgeR << "," << line->UpEdgeR << ") --> " << line->binR << " chs" << endl;
		logfile << " Right compton counts = " << line->countsR << endl;
		logfile << " Estimated Compton counts in signal region = " << line->CompCnts << " +- " << line->CompCnts_err << endl;
		logfile << " Estimated net line counts = " << line->LineCnts << " +- " << line->LineCnts_err << endl;
		logfile << " Sample detection limit = " << line->LdCnts << endl;
		if(line->Detected) logfile << " Line detected! " << endl;
		logfile << "------------" << endl << endl << endl;
	}
	
	logfile << "Measure life time: " << md_aqtime << " s = " << md_aqtime/(24.*3600) << " d" << endl;
	logfile << "Background life time: " << md_Bgtime << " s = " << md_Bgtime/(24.*3600) << " d" << endl;
	
	logfile.close();
	
	mb_dataAnalysed = true;
	
	return;
	
}


void GatorSampleAnalysis::WriteEffTable(string txtfilename){
	
	if(!mb_linesloaded){
		cerr << "GatorSampleAnalysis::WriteEffTable --> ERROR::: No lines loaded!" << endl;
		return;
	}
	
	int nLines = GetNLines();
	
	ofstream efftabfile(txtfilename.c_str());
	
	for(int iLine=0; iLine<nLines; iLine++){
		LineStruct *line = GetLine(iLine);
		
		efftabfile << "|  <sup>" << line->mass << "</sup>" << line->element << "  |  " << line->Energy << "  |  " << line->BR << "  |  ";
		efftabfile << formatdigits2(line->Effic) << "  |  ";
		efftabfile << formatdigits2(line->BRxEffic) << "  |  |" << endl;
		
	}
	
	efftabfile.close();
	
	return;
}


void GatorSampleAnalysis::WriteOutputTable(string outfilename){
	
	if(!mb_dataAnalysed){
		doAnalysis();
		if(!mb_dataAnalysed){
			cerr << "GatorSampleAnalysis::WriteOutputTable --> Data not analysed. Cannot produce any output table." << endl;
			return;
		}
	}
	
	int nLines = GetNLines();
	
	ofstream outfile(outfilename.c_str());
	
	string tmpstr;
	
	for(int iLine=0; iLine<nLines; iLine++){
		LineStruct *line = GetLine(iLine);
		
		outfile << "|  <sup>" << line->mass << "</sup>" << line->element << "  |  " << line->Energy << "  |  ";
		
		tmpstr = formatdigits1(line->PeakCnts,line->PeakCnts_err);
		outfile << tmpstr << "  |  ";
		
		tmpstr = formatdigits1(line->CompCnts,line->CompCnts_err);
		outfile << tmpstr << "  |  ";
		
		if(mb_BGflag){
			tmpstr = formatdigits1((md_aqtime/md_Bgtime)*line->BgLineCnts,(md_aqtime/md_Bgtime)*line->BgLineCnts_err);
			outfile << tmpstr << "  |  ";
		}
		
		tmpstr = formatdigits1(line->LineCnts,line->LineCnts_err);
		outfile << tmpstr << "  |  ";
		
		tmpstr = formatdigits2(line->LdCnts);
		outfile << tmpstr << "  |  ";
		
		tmpstr = formatdigits2(line->LdActiv);
		outfile << tmpstr << "  |  ";
		
		if(line->Detected){
			tmpstr = formatdigits1(line->Activ,line->Activ_err);
			outfile << tmpstr << "  |" << endl;
		}else{
			if(line->Activ<=0){
				tmpstr = formatdigits2(line->LdActiv);
				outfile << "< " << tmpstr << "  |" << endl;
			}else{
				tmpstr = formatdigits2(line->LdActiv + line->Activ);
				outfile << "< " << tmpstr << "  |" << endl;
			}
		}
	}
	
	outfile.close();
	
}


void GatorSampleAnalysis::WriteRootFile(string filename, string mode){
	
	bool recreate_flag,update_flag;
	
	if( (mode!=string("update")) && (mode!=string("UPDATE")) ){
		if ( (mode!=string("recreate")) && (mode!=string("RECREATE")) ){
			cerr << "GatorSampleAnalysis::WriteRootFile --> The saving mode can only be \"UPDATE\" or \"RECREATE\"." << endl;
			return;
		}else{
			update_flag = false;
			recreate_flag = true;
		}
	}else{
		update_flag = true;
		recreate_flag = false;
	}
	
	if(!mb_dataAnalysed){
		doAnalysis();
		if(!mb_dataAnalysed){
			cerr << "GatorSampleAnalysis::WriteRootFile --> Data not analysed. Cannot produce any output file." << endl;
			return;
		}
	}
	
	
	
	TFile *outfile;
	
	TTree *t1;
	
	LineStruct tmpLine;
	string *pstr_mass = &tmpLine.mass;
	string *pstr_elem = &tmpLine.element;
	
	if( recreate_flag ){
		
		outfile = new TFile(filename.c_str(), "recreate");
		t1 = new TTree("t1","Screening results");
		
		t1 -> Branch("mass","string",&pstr_mass);
		t1 -> Branch("element","string",&pstr_elem);
		t1 -> Branch("energy",&tmpLine.Energy,"energy/D");
		t1 -> Branch("energy_err",&tmpLine.Energy_err,"energy_err/D");
		t1 -> Branch("BR",&tmpLine.BR,"BR/D");
		t1 -> Branch("effic",&tmpLine.Effic,"effic/D");
		t1 -> Branch("BRxEffic",&tmpLine.BRxEffic,"BRxEffic/D");
		
		t1 -> Branch("MCAcenter",&tmpLine.MCAcenter,"MCAcenter/D");
		t1 -> Branch("MCAsigma",&tmpLine.MCAsigma,"MCAsigma/D");
		t1 -> Branch("LowEdgeL",&tmpLine.LowEdgeL,"LowEdgeL/D");
		t1 -> Branch("UpEdgeL",&tmpLine.UpEdgeL,"UpEdgeL/D");
		t1 -> Branch("LowEdgeS",&tmpLine.LowEdgeS,"LowEdgeS/D");
		t1 -> Branch("UpEdgeS",&tmpLine.UpEdgeS,"UpEdgeS/D");
		t1 -> Branch("LowEdgeR",&tmpLine.LowEdgeR,"LowEdgeR/D");
		t1 -> Branch("UpEdgeR",&tmpLine.UpEdgeR,"UpEdgeR/D");
		t1 -> Branch("binL",&tmpLine.binL,"binL/D");
		t1 -> Branch("countsL",&tmpLine.countsL,"countsL/D");
		t1 -> Branch("binS",&tmpLine.binS,"binS/D");
		t1 -> Branch("countsS",&tmpLine.countsS,"countsS/D");
		t1 -> Branch("binR",&tmpLine.binR,"binR/D");
		t1 -> Branch("countsR",&tmpLine.countsR,"countsR/D");
		t1 -> Branch("PeakCnts",&tmpLine.PeakCnts,"PeakCnts/D");
		t1 -> Branch("PeakCnts_err",&tmpLine.PeakCnts_err,"PeakCnts_err/D");
		t1 -> Branch("CompCnts",&tmpLine.CompCnts,"CompCnts/D");
		t1 -> Branch("CompCnts_err",&tmpLine.CompCnts_err,"CompCnts_err/D");
		t1 -> Branch("LineCnts",&tmpLine.LineCnts,"LineCnts/D");
		t1 -> Branch("LineCnts_err",&tmpLine.LineCnts_err,"LineCnts_err/D");
		t1 -> Branch("LdCnts",&tmpLine.LdCnts,"LdCnts/D");
		t1 -> Branch("Detected",&tmpLine.Detected,"Detected/O");
		
		if(mb_BGflag){
			t1 -> Branch("BgMCAcenter",&tmpLine.BgMCAcenter,"BgMCAcenter/D");
			t1 -> Branch("BgMCAsigma",&tmpLine.BgMCAsigma,"BgMCAsigma/D");
			t1 -> Branch("BgLowEdgeL",&tmpLine.BgLowEdgeL,"BgLowEdgeL/D");
			t1 -> Branch("BgUpEdgeL",&tmpLine.BgUpEdgeL,"UpBgEdgeL/D");
			t1 -> Branch("BgLowEdgeS",&tmpLine.BgLowEdgeS,"BgLowEdgeS/D");
			t1 -> Branch("BgUpEdgeS",&tmpLine.BgUpEdgeS,"BgUpEdgeS/D");
			t1 -> Branch("BgLowEdgeR",&tmpLine.BgLowEdgeR,"BgLowEdgeR/D");
			t1 -> Branch("BgUpEdgeR",&tmpLine.BgUpEdgeR,"BgUpEdgeR/D");
			t1 -> Branch("BgbinL",&tmpLine.BgbinL,"BgbinL/D");
			t1 -> Branch("BgcountsL",&tmpLine.BgcountsL,"BgcountsL/D");
			t1 -> Branch("BgbinS",&tmpLine.BgbinS,"BgbinS/D");
			t1 -> Branch("BgcountsS",&tmpLine.BgcountsS,"BgcountsS/D");
			t1 -> Branch("BgbinR",&tmpLine.BgbinR,"BgbinR/D");
			t1 -> Branch("BgcountsR",&tmpLine.BgcountsR,"BgcountsR/D");
			t1 -> Branch("BgPeakCnts",&tmpLine.BgPeakCnts,"BgPeakCnts/D");
			t1 -> Branch("BgPeakCnts_err",&tmpLine.BgPeakCnts_err,"BgPeakCnts_err/D");
			t1 -> Branch("BgCompCnts",&tmpLine.BgCompCnts,"BgCompCnts/D");
			t1 -> Branch("BgCompCnts_err",&tmpLine.BgCompCnts_err,"BgCompCnts_err/D");
			t1 -> Branch("BgLineCnts",&tmpLine.BgLineCnts,"BgLineCnts/D");
			t1 -> Branch("BgLineCnts_err",&tmpLine.BgLineCnts_err,"BgLineCnts_err/D");
			t1 -> Branch("BgLdCnts",&tmpLine.BgLdCnts,"BgLdCnts/D");
			t1 -> Branch("BgDetected",&tmpLine.BgDetected,"BgDetected/O");
		}
		
		t1 -> Branch("LdActiv",&tmpLine.LdActiv,"LdActiv/D");
		t1 -> Branch("Activ",&tmpLine.Activ,"Activ/D");
		t1 -> Branch("Activ_err",&tmpLine.Activ_err,"Activ_err/D");
		
		int nLines = GetNLines();
		
		for(int iLine=0; iLine<nLines; iLine++){
			tmpLine = *(GetLine(iLine));
			t1 -> Fill();
		}
		
		outfile->WriteTObject(t1, 0, "overwrite");
		outfile->Close();
		if(outfile) delete outfile;
		if(t1) delete t1;
		
	}else if(update_flag){
		vector<LineStruct> v_tmpLines;
		outfile = new TFile(filename.c_str(), "read");
		t1 = (TTree*)outfile->Get("t1");
		
		t1 -> SetBranchAddress("mass",&pstr_mass);
		t1 -> SetBranchAddress("element",&pstr_elem);
		t1 -> SetBranchAddress("energy",&tmpLine.Energy);
		t1 -> SetBranchAddress("energy_err",&tmpLine.Energy_err);
		t1 -> SetBranchAddress("BR",&tmpLine.BR);
		t1 -> SetBranchAddress("effic",&tmpLine.Effic);
		t1 -> SetBranchAddress("BRxEffic",&tmpLine.BRxEffic);
		
		t1 -> SetBranchAddress("MCAcenter",&tmpLine.MCAcenter);
		t1 -> SetBranchAddress("MCAsigma",&tmpLine.MCAsigma);
		t1 -> SetBranchAddress("LowEdgeL",&tmpLine.LowEdgeL);
		t1 -> SetBranchAddress("UpEdgeL",&tmpLine.UpEdgeL);
		t1 -> SetBranchAddress("LowEdgeS",&tmpLine.LowEdgeS);
		t1 -> SetBranchAddress("UpEdgeS",&tmpLine.UpEdgeS);
		t1 -> SetBranchAddress("LowEdgeR",&tmpLine.LowEdgeR);
		t1 -> SetBranchAddress("UpEdgeR",&tmpLine.UpEdgeR);
		t1 -> SetBranchAddress("binL",&tmpLine.binL);
		t1 -> SetBranchAddress("countsL",&tmpLine.countsL);
		t1 -> SetBranchAddress("binS",&tmpLine.binS);
		t1 -> SetBranchAddress("countsS",&tmpLine.countsS);
		t1 -> SetBranchAddress("binR",&tmpLine.binR);
		t1 -> SetBranchAddress("countsR",&tmpLine.countsR);
		t1 -> SetBranchAddress("PeakCnts",&tmpLine.PeakCnts);
		t1 -> SetBranchAddress("PeakCnts_err",&tmpLine.PeakCnts_err);
		t1 -> SetBranchAddress("CompCnts",&tmpLine.CompCnts);
		t1 -> SetBranchAddress("CompCnts_err",&tmpLine.CompCnts_err);
		t1 -> SetBranchAddress("LineCnts",&tmpLine.LineCnts);
		t1 -> SetBranchAddress("LineCnts_err",&tmpLine.LineCnts_err);
		t1 -> SetBranchAddress("LdCnts",&tmpLine.LdCnts);
		t1 -> SetBranchAddress("Detected",&tmpLine.Detected);
		
		if(mb_BGflag){
			t1 -> SetBranchAddress("BgMCAcenter",&tmpLine.BgMCAcenter);
			t1 -> SetBranchAddress("BgMCAsigma",&tmpLine.BgMCAsigma);
			t1 -> SetBranchAddress("BgLowEdgeL",&tmpLine.BgLowEdgeL);
			t1 -> SetBranchAddress("BgUpEdgeL",&tmpLine.BgUpEdgeL);
			t1 -> SetBranchAddress("BgLowEdgeS",&tmpLine.BgLowEdgeS);
			t1 -> SetBranchAddress("BgUpEdgeS",&tmpLine.BgUpEdgeS);
			t1 -> SetBranchAddress("BgLowEdgeR",&tmpLine.BgLowEdgeR);
			t1 -> SetBranchAddress("BgUpEdgeR",&tmpLine.BgUpEdgeR);
			t1 -> SetBranchAddress("BgbinL",&tmpLine.BgbinL);
			t1 -> SetBranchAddress("BgcountsL",&tmpLine.BgcountsL);
			t1 -> SetBranchAddress("BgbinS",&tmpLine.BgbinS);
			t1 -> SetBranchAddress("BgcountsS",&tmpLine.BgcountsS);
			t1 -> SetBranchAddress("BgbinR",&tmpLine.BgbinR);
			t1 -> SetBranchAddress("BgcountsR",&tmpLine.BgcountsR);
			t1 -> SetBranchAddress("BgPeakCnts",&tmpLine.BgPeakCnts);
			t1 -> SetBranchAddress("BgPeakCnts_err",&tmpLine.BgPeakCnts_err);
			t1 -> SetBranchAddress("BgCompCnts",&tmpLine.BgCompCnts);
			t1 -> SetBranchAddress("BgCompCnts_err",&tmpLine.BgCompCnts_err);
			t1 -> SetBranchAddress("BgLineCnts",&tmpLine.BgLineCnts);
			t1 -> SetBranchAddress("BgLineCnts_err",&tmpLine.BgLineCnts_err);
			t1 -> SetBranchAddress("BgLdCnts",&tmpLine.BgLdCnts);
			t1 -> SetBranchAddress("BgDetected",&tmpLine.BgDetected);
		}
		
		t1 -> SetBranchAddress("LdActiv",&tmpLine.LdActiv);
		t1 -> SetBranchAddress("Activ",&tmpLine.Activ);
		t1 -> SetBranchAddress("Activ_err",&tmpLine.Activ_err);
		
		//Fill the vector with the lines already present in the file
		int nEntries = t1->GetEntries();
		for(int iEnt=0; iEnt<nEntries; iEnt++){
			t1 -> GetEntry(iEnt);
			v_tmpLines.push_back(tmpLine);
		}
		
		outfile->Close();
		if(outfile) delete outfile;
		if(t1) delete t1;
		
		//Extend the vector with the new lines
		int nLines = GetNLines();
		for(int iLine=0; iLine<nLines; iLine++){
			tmpLine = *(GetLine(iLine));
			v_tmpLines.push_back(tmpLine);
		}
		
		//Update the total number of lines
		nLines = v_tmpLines.size();
		
		//Recreate the file and a new tree
		outfile = new TFile(filename.c_str(), "recreate");
		t1 = new TTree("t1","Screening results");
		
		t1 -> Branch("mass","string",&pstr_mass);
		t1 -> Branch("element","string",&pstr_elem);
		t1 -> Branch("energy",&tmpLine.Energy,"energy/D");
		t1 -> Branch("energy_err",&tmpLine.Energy_err,"energy_err/D");
		t1 -> Branch("BR",&tmpLine.BR,"BR/D");
		t1 -> Branch("effic",&tmpLine.Effic,"effic/D");
		t1 -> Branch("BRxEffic",&tmpLine.BRxEffic,"BR_effic/D");
		
		t1 -> Branch("MCAcenter",&tmpLine.MCAcenter,"MCAcenter/D");
		t1 -> Branch("MCAsigma",&tmpLine.MCAsigma,"MCAsigma/D");
		t1 -> Branch("LowEdgeL",&tmpLine.LowEdgeL,"LowEdgeL/D");
		t1 -> Branch("UpEdgeL",&tmpLine.UpEdgeL,"UpEdgeL/D");
		t1 -> Branch("LowEdgeS",&tmpLine.LowEdgeS,"LowEdgeS/D");
		t1 -> Branch("UpEdgeS",&tmpLine.UpEdgeS,"UpEdgeS/D");
		t1 -> Branch("LowEdgeR",&tmpLine.LowEdgeR,"LowEdgeR/D");
		t1 -> Branch("UpEdgeR",&tmpLine.UpEdgeR,"UpEdgeR/D");
		t1 -> Branch("binL",&tmpLine.binL,"binL/D");
		t1 -> Branch("countsL",&tmpLine.countsL,"countsL/D");
		t1 -> Branch("binS",&tmpLine.binS,"binS/D");
		t1 -> Branch("countsS",&tmpLine.countsS,"countsS/D");
		t1 -> Branch("binR",&tmpLine.binR,"binR/D");
		t1 -> Branch("countsR",&tmpLine.countsR,"countsR/D");
		t1 -> Branch("PeakCnts",&tmpLine.PeakCnts,"PeakCnts/D");
		t1 -> Branch("PeakCnts_err",&tmpLine.PeakCnts_err,"PeakCnts_err/D");
		t1 -> Branch("CompCnts",&tmpLine.CompCnts,"CompCnts/D");
		t1 -> Branch("CompCnts_err",&tmpLine.CompCnts_err,"CompCnts_err/D");
		t1 -> Branch("LineCnts",&tmpLine.LineCnts,"LineCnts/D");
		t1 -> Branch("LineCnts_err",&tmpLine.LineCnts_err,"LineCnts_err/D");
		t1 -> Branch("LdCnts",&tmpLine.LdCnts,"LdCnts/D");
		t1 -> Branch("Detected",&tmpLine.Detected,"Detected/O");
		
		if(mb_BGflag){
			t1 -> Branch("BgMCAcenter",&tmpLine.BgMCAcenter,"BgMCAcenter/D");
			t1 -> Branch("BgMCAsigma",&tmpLine.BgMCAsigma,"BgMCAsigma/D");
			t1 -> Branch("BgLowEdgeL",&tmpLine.BgLowEdgeL,"BgLowEdgeL/D");
			t1 -> Branch("BgUpEdgeL",&tmpLine.BgUpEdgeL,"UpBgEdgeL/D");
			t1 -> Branch("BgLowEdgeS",&tmpLine.BgLowEdgeS,"BgLowEdgeS/D");
			t1 -> Branch("BgUpEdgeS",&tmpLine.BgUpEdgeS,"BgUpEdgeS/D");
			t1 -> Branch("BgLowEdgeR",&tmpLine.BgLowEdgeR,"BgLowEdgeR/D");
			t1 -> Branch("BgUpEdgeR",&tmpLine.BgUpEdgeR,"BgUpEdgeR/D");
			t1 -> Branch("BgbinL",&tmpLine.BgbinL,"BgbinL/D");
			t1 -> Branch("BgcountsL",&tmpLine.BgcountsL,"BgcountsL/D");
			t1 -> Branch("BgbinS",&tmpLine.BgbinS,"BgbinS/D");
			t1 -> Branch("BgcountsS",&tmpLine.BgcountsS,"BgcountsS/D");
			t1 -> Branch("BgbinR",&tmpLine.BgbinR,"BgbinR/D");
			t1 -> Branch("BgcountsR",&tmpLine.BgcountsR,"BgcountsR/D");
			t1 -> Branch("BgPeakCnts",&tmpLine.BgPeakCnts,"BgPeakCnts/D");
			t1 -> Branch("BgPeakCnts_err",&tmpLine.BgPeakCnts_err,"BgPeakCnts_err/D");
			t1 -> Branch("BgCompCnts",&tmpLine.BgCompCnts,"BgCompCnts/D");
			t1 -> Branch("BgCompCnts_err",&tmpLine.BgCompCnts_err,"BgCompCnts_err/D");
			t1 -> Branch("BgLineCnts",&tmpLine.BgLineCnts,"BgLineCnts/D");
			t1 -> Branch("BgLineCnts_err",&tmpLine.BgLineCnts_err,"BgLineCnts_err/D");
			t1 -> Branch("BgLdCnts",&tmpLine.BgLdCnts,"BgLdCnts/D");
			t1 -> Branch("BgDetected",&tmpLine.BgDetected,"BgDetected/O");
		}
		
		t1 -> Branch("LdActiv",&tmpLine.LdActiv,"LdActiv/D");
		t1 -> Branch("Activ",&tmpLine.Activ,"Activ/D");
		t1 -> Branch("Activ_err",&tmpLine.Activ_err,"Activ_err/D");
		
		//Fill the new tree
		for(int iLine=0; iLine<nLines; iLine++){
			tmpLine = v_tmpLines.at(iLine);
			t1 -> Fill();
		}
		
		outfile->WriteTObject(t1, 0, "overwrite");
		outfile->Close();
		if(outfile) delete outfile;
		if(t1) delete t1;
	}
	
	
	
}


void GatorSampleAnalysis::DrawActivityPlots(TPad* pad){
	
	if(pad) pad->cd();
	pad -> SetLogy();
	
	if(!mb_dataAnalysed){
		doAnalysis();
		if(!mb_dataAnalysed){
			cerr << "GatorSampleAnalysis::MakeActivityPlots --> Data not analysed. Cannot produce any activity plot." << endl;
			return;
		}
	}
	
	gr_DetActiv = new TGraphErrors();
	gr_DetActiv -> SetMarkerStyle(20);
	gr_DetActiv -> SetMarkerColor(kRed);
	gr_DetActiv -> SetLineColor(kRed);
	
	gr_ULActiv = new TGraphErrors();
	gr_ULActiv -> SetMarkerStyle(20);
	gr_ULActiv -> SetMarkerColor(kBlue);
	
	
	int nLines = GetNLines();
	
	int act_counter = 0;
	int lim_counter = 0;
	for(int iLine=0; iLine<nLines; iLine++){
		LineStruct *line = GetLine(iLine);
		
		if(line->Detected){
			gr_DetActiv -> SetPoint(act_counter,line->Energy,line->Activ);
			gr_DetActiv -> SetPointError(act_counter,0,line->Activ_err);
			act_counter++;
		}else{
			if(line->Activ<=0){
				gr_ULActiv -> SetPoint(lim_counter,line->Energy,line->LdActiv);
				lim_counter++;
			}else{
				gr_ULActiv -> SetPoint(lim_counter,line->Energy,line->LdActiv+line->Activ);
				lim_counter++;
			}
		}
		
	}
	
	TH2F* graphsupport;
	
	if(gr_DetActiv->GetN()==0){
		graphsupport = makeGraphSupport(0,0,TMath::MinElement(gr_ULActiv->GetN(),gr_ULActiv->GetY()),TMath::MaxElement(gr_ULActiv->GetN(),gr_ULActiv->GetY()));
	}else if(gr_ULActiv->GetN()==0){
		graphsupport = makeGraphSupport(TMath::MinElement(gr_DetActiv->GetN(),gr_DetActiv->GetY()),TMath::MaxElement(gr_DetActiv->GetN(),gr_DetActiv->GetY()),0,0);
	}else{
		graphsupport = makeGraphSupport( TMath::MinElement(gr_DetActiv->GetN(),gr_DetActiv->GetY()), TMath::MaxElement(gr_DetActiv->GetN(),gr_DetActiv->GetY()), TMath::MinElement(gr_ULActiv->GetN(),gr_ULActiv->GetY()), TMath::MaxElement(gr_ULActiv->GetN(),gr_ULActiv->GetY()) );
	}
	
	graphsupport -> GetXaxis() -> SetTitle( "Energy [keV]" );
	graphsupport -> GetYaxis() -> SetTitle( (string("mBq/")+m_unit).c_str() );
	graphsupport -> Draw();
	
	if(gr_DetActiv->GetN()>0)gr_DetActiv->Draw("P");
	if(gr_ULActiv->GetN()>0)gr_ULActiv->Draw("P");
	
	TLegend* leg2 = new TLegend(0.60,0.90,0.90,0.80);
	leg2 -> AddEntry(gr_DetActiv,"Detected","P");
	leg2 -> AddEntry(gr_ULActiv,"Upper Limits","P");
	leg2 -> Draw();
}







TH2F* makeGraphSupport(double fMinimum1,double fMaximum1,double fMinimum2,double fMaximum2){
	
	Double_t maxvalue=0;
	Double_t minvalue=0;
	
	/*
	for(Int_t i=0; i<v_vals.size(); i++){
		
		if(i==0){
			maxvalue = v_vals[i];
			minvalue = maxvalue;
		}else{
			if(v_vals[i] > maxvalue) maxvalue = v_vals[i];
			if(v_vals[i] < minvalue) minvalue = v_vals[i];
		}
	}*/
	
	if(fMaximum1>fMaximum2){
		
		maxvalue = fMaximum1;
	}else{
		maxvalue = fMaximum2;
	}
	
	if(fMinimum1<fMinimum2){
		minvalue = fMinimum1;
	}else{
		minvalue = fMinimum2;
	}
	
	//cout << "minvalue = " << minvalue << endl;
	//cout << "maxvalue = " << maxvalue << endl;
	
	TH2F* graphsupport = new TH2F("graphsupport","",1000,0.,2700.,1000,0.,maxvalue+0.1*(maxvalue-minvalue));
	
	return graphsupport;
	
} 
