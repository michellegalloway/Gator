#ifndef GATOR_CALIB_CLASS_HH
#define GATOR_CALIB_CLASS_HH

/*
The classes in this header are designed to be used interactively with root cling interpreter
*/

#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"


using namespace std;

namespace Gator{
	
#include "GatorStructs.h"
	
	class GatorCalib{
	public:
		GatorCalib();
		
		void LoadCalibFiles(const string& sourcename, const string& dir);
		TH1D* SelectSpectrum(const string& sourcename);
		TH1D* DrawSpectrum(const string& opt=string(""));
		void AddLine(const CalibLine *line);
		
		void LoadLinesFromTree(const string& rootfile);
		
		TH1D* GetSpectrum(const string& sourcename);
		CalibLine* GetCalibLine(const string& linename);
		
		BCHistogramFitter* FitLine(const string& linename);
		
		bool SaveLines(const string& rootfile);
		
	private:
		
		double peakFitFunc(double* x, double* par);
		
		void amplInit(TH1D* histo, CalibLine& line);
		void costInit(TH1D* histo, CalibLine& line);
		void stepInit(TH1D* histo, CalibLine& line);
		void sigmaInit(TH1D* histo, CalibLine& line);
		
		map<string, CalibLine*> fLinesmap;
		
		map<string, TH1D*> fSpectramap;
		
		TH1D* fHisto;
	};
	
}

#endif