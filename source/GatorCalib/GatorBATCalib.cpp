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

#include "GatorStructs.h"
#include "GatorCalibFunc.h"
#include "GatorCalibScr.h"
#include "misc.h"

using namespace std;

void usage();

extern TApplication* theApp=NULL;

int main(int argc, char* argv[]){
	
	int mickey = 0;
	char *mouse[10];
	theApp = new TApplication("theApp", &mickey, mouse );
	
	char c;
	
	bool calibset_flag = false;
	bool configfile_flag = false;
	bool recreate = false; //To recreat a file
	bool update = false; //To force the update mode
	
	string calibset("");
	string configfile("");
	
    // parse switches
	while((c = getopt(argc,argv,"d:f:ruh")) != -1)
	{
		switch(c)
		{
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
			update = false;
			break;
			
			case 'u':
			update = true;
			recreate = false;
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
	
	return GatorCalibScriptBAT(calibset, configfile, recreate, update);
}

void usage(){return;}
