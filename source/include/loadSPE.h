#ifndef LOADSPE_H
#define LOADSPE_H

#include "TH1D.h"

#include <string>


TH1D* loadSpe(const char* dir, Double_t& aqtime);//This is to load the calibration SPE files

TH1D* loadSPE(const char* dir, Double_t& aqtime);//This is to load the SPE files from standard acquisition

TH1D* loadSPE3(Char_t* filename, Double_t& time, UInt_t unixtime, std::string* descr_str);

TH1D* loadSingleSPE(const char* name, Double_t& aqtime);


#endif
