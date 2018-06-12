#include <cstdlib>

#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "TH1D.h"

#include "GatorCounter.hh"
#include "GatorStructs.h"


using namespace std;

//Default Constructor
GatorCounter::GatorCounter(){
	
	mb_BGflag = true;
	
	mb_nonIntCounts = true;
	
	mp_sampleSpect = NULL;
	mp_backgroundSpect = NULL;
	
	mv_lines = new vector<LineStruct>;
	
}


int GatorCounter::GetNLines(){
	
	return mv_lines->size();
	
}


LineStruct* GatorCounter::GetLine(int iLine){
	
	if( (iLine>=mv_lines->size()) || (iLine<0) ) return NULL;
	
	return &(mv_lines->at(iLine));
	
}


void GatorCounter::LineCounts(LineStruct &line){
	
	
	if(mb_nonIntCounts){
		
		NonIntRegionCounts(line);
		
	}else{
		
		PoissonRegionCounts(line);
		
	}
	
	return;
}


void GatorCounter::NonIntRegionCounts(LineStruct &line){
	
	int bin1, bin2, bin3, bin4, bin5, bin6;
	
	//Fraction of counts to be considered in the edge bis
	double fracUp, fracLow, BinLowEdge, BinWidth;
	
	
	//Sample spectrum
	bin1 = mp_sampleSpect -> FindBin(line.LowEdgeL);
	bin2 = mp_sampleSpect -> FindBin(line.UpEdgeL);
	bin3 = mp_sampleSpect -> FindBin(line.LowEdgeS);
	bin4 = mp_sampleSpect -> FindBin(line.UpEdgeS);
	bin5 = mp_sampleSpect -> FindBin(line.LowEdgeR);
	bin6 = mp_sampleSpect -> FindBin(line.UpEdgeR);
	
	/*
	cout << "\nbin1 = " << bin1 << endl;
	cout << "bin2 = " << bin2 << endl;
	cout << "bin3 = " << bin3 << endl;
	cout << "bin4 = " << bin4 << endl;
	cout << "bin5 = " << bin5 << endl;
	cout << "bin6 = " << bin6 << endl;
	exit(0);
	*/
	
	//Left compton region
	line.binL = line.UpEdgeL - line.LowEdgeL;
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin1);
	BinWidth = mp_sampleSpect->GetBinWidth(bin1);
	
	fracLow = 1.-(line.LowEdgeL-BinLowEdge)/BinWidth;
	
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin2);
	BinWidth = mp_sampleSpect->GetBinWidth(bin2);
	
	fracUp = (line.UpEdgeL-BinLowEdge)/BinWidth;
	
	line.countsL = fracLow*mp_sampleSpect->GetBinContent(bin1);
	for(int iBin=bin1+1; iBin<bin2; iBin++){
		line.countsL += mp_sampleSpect->GetBinContent(iBin);
	}
	line.countsL += fracUp*mp_sampleSpect->GetBinContent(bin2);
	
	//Signal region
	line.binS = line.UpEdgeS - line.LowEdgeS;
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin3);
	BinWidth = mp_sampleSpect->GetBinWidth(bin3);
	
	fracLow = 1.-(line.LowEdgeS-BinLowEdge)/BinWidth;
	
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin4);
	BinWidth = mp_sampleSpect->GetBinWidth(bin4);
	
	fracUp = (line.UpEdgeS-BinLowEdge)/BinWidth;
	
	line.countsS = fracLow*mp_sampleSpect->GetBinContent(bin3);
	for(int iBin=bin3+1; iBin<bin4; iBin++){
		line.countsS += mp_sampleSpect->GetBinContent(iBin);
	}
	line.countsS += fracUp*mp_sampleSpect->GetBinContent(bin4);
	
	//Right compton region
	line.binR = line.UpEdgeR - line.LowEdgeR;
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin5);
	BinWidth = mp_sampleSpect->GetBinWidth(bin5);
	
	fracLow = 1.-(line.LowEdgeR-BinLowEdge)/BinWidth;
	
	BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin6);
	BinWidth = mp_sampleSpect->GetBinWidth(bin6);
	
	fracUp = (line.UpEdgeR-BinLowEdge)/BinWidth;
	
	line.countsR = fracLow*mp_sampleSpect->GetBinContent(bin5);
	for(int iBin=bin5+1; iBin<bin6; iBin++){
		line.countsR += mp_sampleSpect->GetBinContent(iBin);
	}
	line.countsR += fracUp*mp_sampleSpect->GetBinContent(bin6);
	
	
	
	//Baground spectrum
	if(mb_BGflag){
		
		bin1 = mp_backgroundSpect -> FindBin(line.BgLowEdgeL);
		bin2 = mp_backgroundSpect -> FindBin(line.BgUpEdgeL);
		bin3 = mp_backgroundSpect -> FindBin(line.BgLowEdgeS);
		bin4 = mp_backgroundSpect -> FindBin(line.BgUpEdgeS);
		bin5 = mp_backgroundSpect -> FindBin(line.BgLowEdgeR);
		bin6 = mp_backgroundSpect -> FindBin(line.BgUpEdgeR);
		
		/*
		cout << "\nbin1 = " << bin1 << endl;
		cout << "bin2 = " << bin2 << endl;
		cout << "bin3 = " << bin3 << endl;
		cout << "bin4 = " << bin4 << endl;
		cout << "bin5 = " << bin5 << endl;
		cout << "bin6 = " << bin6 << endl;
		exit(0);
		*/
		
		//Left compton region
		line.BgbinL = line.BgUpEdgeL - line.BgLowEdgeL;
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin1);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin1);
	
		fracLow = 1.-(line.BgLowEdgeL-BinLowEdge)/BinWidth;
	
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin2);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin2);
	
		fracUp = (line.BgUpEdgeL-BinLowEdge)/BinWidth;
	
		line.BgcountsL = fracLow*mp_backgroundSpect->GetBinContent(bin1);
		for(int iBin=bin1+1; iBin<bin2; iBin++){
			line.BgcountsL += mp_backgroundSpect->GetBinContent(iBin);
		}
		line.BgcountsL += fracUp*mp_backgroundSpect->GetBinContent(bin2);
	
		//Signal region
		line.BgbinS = line.BgUpEdgeS - line.BgLowEdgeS;
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin3);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin3);
	
		fracLow = 1.-(line.BgLowEdgeS-BinLowEdge)/BinWidth;
	
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin4);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin4);
	
		fracUp = (line.BgUpEdgeS-BinLowEdge)/BinWidth;
	
		line.BgcountsS = fracLow*mp_backgroundSpect->GetBinContent(bin3);
		for(int iBin=bin3+1; iBin<bin4; iBin++){
			line.BgcountsS += mp_backgroundSpect->GetBinContent(iBin);
		}
		line.BgcountsS += fracUp*mp_backgroundSpect->GetBinContent(bin4);
	
		//Right compton region
		line.BgbinR = line.BgUpEdgeR - line.BgLowEdgeR;
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin5);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin5);
	
		fracLow = 1.-(line.BgLowEdgeR-BinLowEdge)/BinWidth;
	
		BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin6);
		BinWidth = mp_backgroundSpect->GetBinWidth(bin6);
	
		fracUp = (line.BgUpEdgeR-BinLowEdge)/BinWidth;
	
		line.BgcountsR = fracLow*mp_backgroundSpect->GetBinContent(bin5);
		for(int iBin=bin5+1; iBin<bin6; iBin++){
			line.BgcountsR += mp_backgroundSpect->GetBinContent(iBin);
		}
		line.BgcountsR += fracUp*mp_backgroundSpect->GetBinContent(bin6);
	}
	
	
}

void GatorCounter::PoissonRegionCounts(LineStruct &line){
	
	int bin1, bin2, bin3, bin4, bin5, bin6;
	
	//Fraction of counts to be considered in the edge bis
	double fracUp, fracLow, BinLowEdge, BinWidth;
	
	
	//Sample spectrum
	bin1 = mp_sampleSpect -> GetXaxis() -> FindBin(line.LowEdgeL);
	bin2 = mp_sampleSpect -> GetXaxis() -> FindBin(line.UpEdgeL);
	bin3 = mp_sampleSpect -> GetXaxis() -> FindBin(line.LowEdgeS);
	bin4 = mp_sampleSpect -> GetXaxis() -> FindBin(line.UpEdgeS);
	bin5 = mp_sampleSpect -> GetXaxis() -> FindBin(line.LowEdgeR);
	bin6 = mp_sampleSpect -> GetXaxis() -> FindBin(line.UpEdgeR);
	
	
	if(bin2==bin3){
		
		BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin2);
		BinWidth = mp_sampleSpect->GetBinWidth(bin2);
		fracUp = (line.UpEdgeL-BinLowEdge)/BinWidth;
		
		BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin3);
		BinWidth = mp_sampleSpect->GetBinWidth(bin3);
		fracLow = 1.-(line.LowEdgeS-BinLowEdge)/BinWidth;
		
		if(fracUp>fracLow){
			bin3 +=1;
		}else{
			bin2 -=1;
		}
	}
	
	if(bin4==bin5){
		
		BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin4);
		BinWidth = mp_sampleSpect->GetBinWidth(bin4);
		fracUp = (line.UpEdgeS-BinLowEdge)/BinWidth;
		
		BinLowEdge = mp_sampleSpect->GetBinLowEdge(bin5);
		BinWidth = mp_sampleSpect->GetBinWidth(bin5);
		fracLow = 1.-(line.LowEdgeR-BinLowEdge)/BinWidth;
		
		if(fracUp>fracLow){
			bin5 +=1;
		}else{
			bin4 -=1;
		}
	}
	
	//Left compton region
	line.binL = (double)(bin2 - bin1 + 1);
	for(int iBin=bin1; iBin<=bin2; iBin++){
		line.countsL += mp_sampleSpect->GetBinContent(iBin);
	}
	
	//Signal region
	line.binS = (double)(bin4 - bin3 + 1);
	for(int iBin=bin3; iBin<=bin4; iBin++){
		line.countsS += mp_sampleSpect->GetBinContent(iBin);
	}
	
	//Right compton region
	line.binR = (double)(bin6 - bin5 + 1);
	for(int iBin=bin5; iBin<=bin6; iBin++){
		line.countsR += mp_sampleSpect->GetBinContent(iBin);
	}
	
	//Baground spectrum
	if(mb_BGflag){
		
		//Sample spectrum
		bin1 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgLowEdgeL);
		bin2 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgUpEdgeL);
		bin3 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgLowEdgeS);
		bin4 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgUpEdgeS);
		bin5 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgLowEdgeR);
		bin6 = mp_backgroundSpect -> GetXaxis() -> FindBin(line.BgUpEdgeR);
	
	
		if(bin2==bin3){
		
			BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin2);
			BinWidth = mp_backgroundSpect->GetBinWidth(bin2);
			fracUp = (line.BgUpEdgeL-BinLowEdge)/BinWidth;
		
			BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin3);
			BinWidth = mp_backgroundSpect->GetBinWidth(bin3);
			fracLow = 1.-(line.BgLowEdgeS-BinLowEdge)/BinWidth;
		
			if(fracUp>fracLow){
				bin3 +=1;
			}else{
				bin2 -=1;
			}
		}
	
		if(bin4==bin5){
		
			BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin4);
			BinWidth = mp_backgroundSpect->GetBinWidth(bin4);
			fracUp = (line.BgUpEdgeS-BinLowEdge)/BinWidth;
		
			BinLowEdge = mp_backgroundSpect->GetBinLowEdge(bin5);
			BinWidth = mp_backgroundSpect->GetBinWidth(bin5);
			fracLow = 1.-(line.BgLowEdgeR-BinLowEdge)/BinWidth;
		
			if(fracUp>fracLow){
				bin5 +=1;
			}else{
				bin4 -=1;
			}
		}
	
		//Left compton region
		line.BgbinL = (double)(bin2 - bin1 + 1);
		for(int iBin=bin1; iBin<=bin2; iBin++){
			line.BgcountsL += mp_backgroundSpect->GetBinContent(iBin);
		}
	
		//Signal region
		line.BgbinS = (double)(bin4 - bin3 + 1);
		for(int iBin=bin3; iBin<=bin4; iBin++){
			line.BgcountsS += mp_backgroundSpect->GetBinContent(iBin);
		}
	
		//Right compton region
		line.BgbinR = (double)(bin6 - bin5 + 1);
		for(int iBin=bin5; iBin<=bin6; iBin++){
			line.BgcountsR += mp_backgroundSpect->GetBinContent(iBin);
		}
		
	}
}