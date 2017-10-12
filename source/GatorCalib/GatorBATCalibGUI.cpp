#include "GatorCalibGUIclass.hh"
#include "GatorCalibClass.hh"
#include "GatorStructs.h"
#include "GatorCalibFunc.h"
#include "GatorCalibScr.h"
#include "misc.h"

#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>

#include <BAT/BCAux.h>
//#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>

#include <string>
#include <cstdlib>
#include <unistd.h>


using namespace std;

void usage();

extern TApplication* theApp=NULL;

int main(int argc, char* argv[]){
	
	theApp = new TApplication("theApp", &argc, argv );
	
	char c;
	
	bool calibset_flag = false;
	bool configfile_flag = false;
	bool recreate = false; //To recreat a file
	bool update = false; //To force the update mode
	bool guimode = false;
	
	string calibset("");
	string configfile("");
	
	// parse switches
	while((c = getopt(argc,argv,"d:f:gruh")) != -1)
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
			
			case 'g':
			guimode = true;
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
	
	
	if(guimode)
	{
		new Gator::GatorCalibGUI;
		theApp.Run();
		return 0;
	}
	
	if(!(calibset_flag && configfile_flag)){
		cout << endl << argv[0] << ":ERROR --> One or more mandatory arguments are missing!!!" << endl << endl;
		usage();
		return(-1);
	}
	
	return GatorCalibScriptBAT(calibset, configfile, recreate, update);
}

void usage(){return;}
