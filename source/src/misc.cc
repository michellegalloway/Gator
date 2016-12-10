#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
//#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>
//#include <TCut.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
//#include <TIterator.h>
//#include <TList.h>
//#include <TMultiGraph.h>
#include <TMath.h>
#include <TApplication.h>
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>
//#include <TLatex.h>
//#include <TTimeStamp.h>
#include <TLine.h>
//#include <TVirtualFitter.h>
//#include <TGaxis.h>
//#include <TMarker.h>
#include <TFitResult.h>


using namespace std;

bool fexist(string filename)
{
	ifstream infile(filename.c_str());
	return infile;
}

string getEnvVar( string const & key )
{
    char * val = getenv( key.c_str() );
    return val == NULL ? string("") : string(val);
}