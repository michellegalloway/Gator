#include "GatorGlobals.hh"
#include "GatorStructs.h"
#include "GatorCalibFunc.h"
#include "GatorCalibScr.h"
#include "misc.h"

#include <string>
#include <cstdlib>
#include <unistd.h>

#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>

#include <BAT/BCAux.h>
//#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>

using namespace std;

void usage();

TApplication* theApp;

int main(int argc, char* argv[]){
	
	int mickey = 0;
	char *mouse[10];
	theApp = new TApplication("theApp", &mickey, mouse );
	
	char c;
	
	bool calibset_flag = false;
	bool configfile_flag = false;
	bool recreate = false;
	
	string calibset("");
	string configfile("");
	
    // parse switches
    while((c = getopt(argc,argv,"d:f:rh")) != -1)
    {
      switch(c)	{
        
        case 'd': 
          calibset_flag = true;
          calibset = optarg;
          break;
        
        case 'f':
          configfile_flag = true;
          configfile = optarg;
          break;
        
        case 'r':
          recreate = true;
          break;
				
        case 'h':
			usage();
			return 0;
          	break;
        
        default:
		cout << endl << argv[0] << ":ERROR --> \"-" << c << " doesn't match with any option!!!" << endl << endl;
          usage();
		  return 0;
      }
    }
	
	if(!(calibset_flag && configfile_flag)){
		cout << endl << argv[0] << ":ERROR --> One or more mandatory arguments are missing!!!" << endl << endl;
		usage();
		return(-1);
	}
	
	return GatorCalibScriptLL(calibset, configfile, recreate);
}

void usage(){return;}
