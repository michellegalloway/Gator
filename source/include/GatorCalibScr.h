#ifndef GATORCALIBSCR_H
#define GATORCALIBSCR_H

#include <string>

int GatorCalibScriptBAT(string calibset, string configfile, bool recreate=false, bool update=false);

int GatorCalibScriptLL(string calibset, string configfile, bool recreate=false);


#endif