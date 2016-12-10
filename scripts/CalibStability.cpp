#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include <TStyle.h>
#include "TFile.h"
//#include <TTree.h>
#include "TH1F.h"
//#include <TH2F.h>
//#include <TH1D.h>
#include <TH2D.h>
#include "TF1.h"
#include "TLegend.h"
//#include <TCut.h>
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include <TGraph.h>
#include <TGraphErrors.h>
//#include <TIterator.h>
//#include <TList.h>
//#include <TMultiGraph.h>
#include "TMath.h"
#include "TApplication.h"
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>
//#include <TLatex.h>
//#include <TTimeStamp.h>
//#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
//#include <TFitResult.h>
#include "TDatime.h"


using namespace std;

class calibPointStruct{
	
public:
	calibPointStruct(string _year, string _month, string _day):
	    year(_year),
	    month(_month),
	    day(_day)
	{
		dataname = year+string("-")+month+string("-")+day;
		
		TDatime dt_tmp( (dataname+string(" 00:00:00")).c_str() );
		unixtime = dt_tmp.Convert();
	};
	~calibPointStruct(){;};
	
	string year;
	string month;
	string day;
	string dataname;
	unsigned unixtime;
	
	
};


void CalibStability(){
	
	gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(10);
    gStyle->SetStatColor(10);
    gStyle->SetStatFont(42);
	gStyle->SetCanvasBorderMode(0);
	
	string calibarchive("/Users/francesco/PhD/Gator/calibrations/archive/");
	
	vector<calibPointStruct> calibpoints;
	calibpoints.push_back( calibPointStruct("2012","06","19") );
	calibpoints.push_back( calibPointStruct("2012","10","04") );
	calibpoints.push_back( calibPointStruct("2013","07","01") );
	calibpoints.push_back( calibPointStruct("2013","12","23") );
	//calibpoints.push_back( calibPointStruct("2014","10","10") );
	
	unsigned nPoints = calibpoints.size();
	
	vector<double> EnergiesVec;
	vector<string> LabelTextVec;
	vector<Color_t> PlotColorVect;
	EnergiesVec.push_back(238.632); LabelTextVec.push_back("239 keV line"); PlotColorVect.push_back(kBlue);
	EnergiesVec.push_back(583.187); LabelTextVec.push_back("583 keV line"); PlotColorVect.push_back(kRed);
	EnergiesVec.push_back(1620.738); LabelTextVec.push_back("1621 keV line"); PlotColorVect.push_back(kGreen+2);
	EnergiesVec.push_back(2614.511); LabelTextVec.push_back("2615 keV line"); PlotColorVect.push_back(kOrange+6);
	
	unsigned nLines = EnergiesVec.size();
	
	vector<TGraph*> PlotsVec(nLines);
	for(unsigned iLine=0; iLine<nLines; iLine++){
		PlotsVec.at(iLine) = new TGraph();
		PlotsVec.at(iLine)->SetMarkerColor(PlotColorVect.at(iLine));
		PlotsVec.at(iLine)->SetMarkerStyle(20);
		PlotsVec.at(iLine)->SetLineColor(PlotColorVect.at(iLine));
		PlotsVec.at(iLine)->SetLineStyle(2);
		PlotsVec.at(iLine)->SetTitle("; Data; MCA channel");
	}
	
	TFile *rootfile = NULL;
	
	double maxMCAch = 0;
	double maxUnixTime, minUnixTime;
	
	for(unsigned iPoint=0; iPoint<nPoints; iPoint++){
		
		string filename = calibarchive+calibpoints.at(iPoint).dataname+string("/calibration.root");
		double timepoint = (double)calibpoints.at(iPoint).unixtime;
		rootfile = TFile::Open(filename.c_str(),"read");
		TF1 *calib_fcn_bw = (TF1*)rootfile->Get("calib_fcn_bw");
		
		if(iPoint==0){
			maxUnixTime = timepoint;
			minUnixTime = timepoint;
		}else{
			if(timepoint>maxUnixTime) maxUnixTime = timepoint;
			if(timepoint<minUnixTime) minUnixTime = timepoint;
		}
		
		for(unsigned iLine=0; iLine<nLines; iLine++){
			double MCAch = calib_fcn_bw->Eval(EnergiesVec.at(iLine));
			PlotsVec.at(iLine)->SetPoint((int)iPoint, timepoint, MCAch);
			if(MCAch>maxMCAch) maxMCAch = MCAch;
		}
		rootfile->Close();
		if(rootfile){
			delete rootfile;
			rootfile=NULL;
		}
	}
	
	
	TCanvas *c1 = new TCanvas("c1","c1",1200,nLines*400);
	c1->Divide(1,nLines);
	
	double Dx = maxUnixTime - minUnixTime;
	TH2D *frameHist1 = NULL;
	
	stringstream ss_tmp;
	
	for(unsigned iLine=0; iLine<nLines; iLine++){
		
		ss_tmp << "frameHist1_" << iLine+1;
		
		int nEntries = PlotsVec.at(iLine)->GetN();
		double *yarr = PlotsVec.at(iLine)->GetY();
		double ymax = TMath::MaxElement( nEntries, yarr );
		double ymin = TMath::MinElement( nEntries, yarr );
		double Dy = ymax-ymin;
		
		frameHist1 = new TH2D( ss_tmp.str().c_str(), "; Data; MCA channel", 200, minUnixTime-Dx/10, maxUnixTime+Dx/10, 200, ymin-Dy/10, ymax+Dy/10);
		
		c1->cd(iLine+1)->SetGridy();
		frameHist1->Draw();
		PlotsVec.at(iLine)->Draw("PL");
		
	}
	
	//Plot in only one frame, all of them normalized by average
	TCanvas *c2 = new TCanvas("c2","c2",1200,400);
	c2->GetGridy();
	
	TLegend *leg = new TLegend(0.60,0.90,0.90,0.80);
	
	double yrelmax, yrelmin;
	
	vector<TGraph*> RelPlotsVec(nLines);
	for(unsigned iLine=0; iLine<nLines; iLine++){
		RelPlotsVec.at(iLine) = new TGraph();
		RelPlotsVec.at(iLine)->SetMarkerColor(PlotColorVect.at(iLine));
		RelPlotsVec.at(iLine)->SetMarkerStyle(20);
		RelPlotsVec.at(iLine)->SetLineColor(PlotColorVect.at(iLine));
		RelPlotsVec.at(iLine)->SetLineStyle(2);
		//RelPlotsVec.at(iLine)->SetTitle("; Data; MCA channel");
		
		int nEntries = PlotsVec.at(iLine)->GetN();
		double *yarr = PlotsVec.at(iLine)->GetY();
		double ymean = TMath::Mean(nEntries, yarr);
		
		for(int iEntry=0; iEntry<nEntries; iEntry++){
			double timepoint = (double)calibpoints.at(iEntry).unixtime;
			double yrel = (yarr[iEntry]/ymean)-1.;
			RelPlotsVec.at(iLine)->SetPoint(iEntry, timepoint, yrel);
			if(iLine==0 && iEntry==0){
				yrelmax = yrel;
				yrelmin = yrel;
			}else{
				if(yrelmax < yrel) yrelmax = yrel;
				if(yrelmin > yrel) yrelmin = yrel;
			}
		}
	}
	
	double Dy = yrelmax-yrelmin;
	
	TH2D *frameHist2 = new TH2D( "frameHist2", "; Data; Relative deviation", 200, minUnixTime-Dx/10, maxUnixTime+Dx/10, 200, yrelmin-Dy/10, yrelmax+Dy/10);
	
	frameHist2->GetXaxis()->SetTimeDisplay(1);
	frameHist2->GetXaxis()->SetTimeFormat("%b\n %Y %F1970-01-01 00:00:00");
	frameHist2->Draw();
	
	for(unsigned iLine=0; iLine<nLines; iLine++){
		leg->AddEntry(RelPlotsVec.at(iLine),LabelTextVec.at(iLine).c_str(),"PL");
		RelPlotsVec.at(iLine)->Draw("PL");
	}
	
	leg -> SetFillColor(kWhite);
	leg -> SetBorderSize(0);
	leg -> SetTextFont(62);
	leg->Draw();
	
	return;
	
}