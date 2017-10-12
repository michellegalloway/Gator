#ifndef GATOR_CALIB_CLASS_HH
#define GATOR_CALIB_CLASS_HH

/*
The classes in this header are designed to be used interactively with root cling interpreter
*/
#include "GatorStructs.h"
#include "loadSPE.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <string>
#include <vector>
#include <map>


using namespace std;

class BCHistogramFitter;

namespace Gator{
	
	class GatorCalib
	{
	protected:
		//Data members for the calibration
		map<string, CalibLine*> fLinesmap;
		map<string, string> fLinesSpectraMap;
		map<string, TH1D*> fSpectramap;
		TH1D* fHisto;
		
		bool fDebug;
		
		static double peakFitFuncA(double* x, double* par);
		static double peakFitFuncB(double* x, double* par);
		
	public:
		GatorCalib();
		virtual ~GatorCalib(){;};
		
		void LoadCalibFiles(const string& sourcename, const string& dir);
		TH1D* SelectSpectrum(const string& sourcename);
		TH1D* DrawSpectrum(const string& opt=string(""));
		void AddLine(CalibLine *line, const string& spectrum);
		void AddLine(const string& _massnum, const string& _element, const double& _litEn, const double& _litEnErr, const string& spectrum);
		void AddLine(const string& _massnum, const string& _element, const double& _litEn, const string& spectrum){
			AddLine(_massnum, _element, _litEn, 0.0, spectrum);
		};
		
		bool LoadLinesFromTree(const string& rootfile);
		
		TH1D* GetSpectrum(const string& sourcename);
		CalibLine* GetCalibLine(const string& linename);
		
		BCHistogramFitter* FitLine(CalibLine& line);
		
		bool SaveLines(const string& rootfile, bool update=true);
		
		void SetDebug(bool flag=true){fDebug=flag;};
		bool IsDebug(){return fDebug;};

	private:
		TH1D* loadSpe(const char* dir, double& aqtime);
		
		
		void amplInit(TH1D* histo, CalibLine& line);
		void costInit(TH1D* histo, CalibLine& line);
		void stepInit(TH1D* histo, CalibLine& line);
		void sigmaInit(TH1D* histo, CalibLine& line);
		
#if defined(__CLING__)
		ClassDef(GatorCalib,0)
#endif
	};

}//End of 'Gator' namespace

#endif
