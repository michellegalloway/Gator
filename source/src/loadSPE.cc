#include "loadSPE.h"

#include "TH1D.h"
#include "TMath.h"
#include "TDatime.h"

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace std;

#if defined(__CLING__)
TH1D* Gator::loadSPE(const char* dir, double& aqtime)
#else
TH1D* loadSPE(const char* dir, double& aqtime)
#endif
{

	//-------------------------------------------------//
	// Load of the sample histogram from the SPE files //
	// this version is for data files                  //
	//-------------------------------------------------//

	string tmpstr("");
	stringstream tmpsstr;
	
	aqtime=0;
	
	Double_t time_tmp=0, entry=0; 
	Int_t channel=0;
	
	vector<Double_t>* entries_tmp = new vector<Double_t>;
	vector<Double_t>* entries = new vector<Double_t>;
	
	string dirname(dir);
	string command("ls ");
	vector<string> spelist;
	
	command += dirname;
	command += string("*.SPE > tmplist");
	
	system(command.c_str());
	ifstream tmplist("tmplist");
	
	while(getline(tmplist,tmpstr)){
		spelist.push_back(tmpstr);
	}
	
	int nSpe = spelist.size();
	
	system("rm -f tmplist");
	 
	ifstream spe;
	
	
	
	for(int iSpe=0; iSpe<nSpe; iSpe++){
		
		string infilename = spelist.at(iSpe);
		cout << "\nLoading file <" << infilename << ">" << endl;
		
		spe.open(infilename.c_str());
		
		if (spe.fail()){
			cout << "Error: can't open the file <" << infilename.c_str() << ">\nSkipping the entry!	" << endl;
			continue;//Skip the loop
		}
		
		entries_tmp->clear();
		
		int channels;
		
		for (Int_t line = 0; line<=11; line++) {
			//Here I go ahead in the file discarding the header lines
			getline(spe,tmpstr);
			
			if (line==9){//Acquisition live time (the first of the two numbers)
				//cout << "Content of \"tmpstr\" at line 9 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				time_tmp = atof(tmpstr.c_str());
				cout << "Acquisition time: " << time_tmp << " secs" << endl;
			}
			
			if (line==11){//Info on how many channels are acquired (the second of the two numbers)
				//cout << "Content of \"tmpstr\" at line 11 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				tmpsstr >> tmpstr;
				channels = atoi(tmpstr.c_str());
				channels++;//The number of channels is 1 more of the last channel number
				cout << "Number of channels: " << channels << endl;
			}
		}
		
		//Loading the ADC counts from the file and put them in the vector entry_tmp
		while (spe.good()) {
			spe >> tmpstr;
			entry = atof(tmpstr.c_str());
			if(iSpe==0){
				entries -> push_back(entry);
				if (entries->size()==channels) break;
			} else {
				entries_tmp -> push_back(entry);
				if (entries_tmp->size()==channels) break;
			}
				
			if (!spe.good()) break;
		}
		
		if(iSpe==0){
			cout << "Data vector size: " << entries->size() << endl;
		}else{
			cout << "Data vector size: " << entries_tmp->size() << endl;
			for(Int_t i=0; i<entries_tmp->size(); i++){
				(*entries)[i] += entries_tmp->at(i);
			}
		}
		
		spe.close();
		
		aqtime += time_tmp; //Stored total total acquisition time
		
	}
	
	//Create and fill the histogram
	TH1D* histoMCA = new TH1D("histoMCA","Spectrum ADC",entries->size(),1,entries_tmp->size());
	for(Int_t i=0; i<entries->size(); i++){
		histoMCA->Fill(i+1,entries->at(i));
	}
	
	//-------------------------------------------------------------//	
	// Finished to load of the sample histogram from the SPE files //
	//-------------------------------------------------------------//

	
	delete entries_tmp;
	delete entries;
	
	return histoMCA;
	
}
//------------------------------------------//


#if defined(__CLING__)
TH1D* Gator::loadSPE3(Char_t* filename, Double_t& time, UInt_t unixtime, string* descr_str)
#else
TH1D* loadSPE3(Char_t* filename, Double_t& time, UInt_t unixtime, string* descr_str)
#endif
{//This pointer can be set to 0 if description is not wanted

	//------------------------------------------------------//	
	// Load of the sample histogram from the SPE files      //
	// This version gives the Unix timestamp and the sample //
	// description if requested                             //
	//------------------------------------------------------//

	string bufferobj;//, dataobj, livetimeobj, channelobj;
	//Char_t buffer[1000]; //Multi-purpose string
	Char_t datastr[20]; //For the date extraction from the spe file (that line is long 19 characters)
	Char_t year[4]="",month[2]="",day[2]="",hour[2]="",min[2]="",sec[2]="";
	
	TDatime td; //Used to put in the acquisition data and to generate the correspondent unix timestamp
	//Int_t timeoffset = -3600; //To be used to return the correct Unix timestamp
	
	ifstream spe(filename);
	
	vector<Double_t>* entries = new vector<Double_t>;
	
	Double_t entry;
	
	Int_t channels=0;
	
	if (spe.fail()){
		cout << "\n\nError: can't open the file <" << filename << ">" << endl;
		return 0;
	}
	
	cout << "\n\nLoading file <" << filename << ">" << endl;
	
	for (Int_t i = 0; i<=11; i++) {//There are 12 header lines
		//Here I go ahead in the file discarding the header lines
		getline(spe,bufferobj);
	
		if (i==1){//In the 2nd line there is the description of the sample
			cout << "Sample description is: " << bufferobj << endl;
			if(!(descr_str==0)){//Take the description only if is requested
				(*descr_str) = bufferobj.substr(0,bufferobj.find("\n"));
			}
		}
		
		if (i==7){//Line with the starting acquisition time to be converted in unix timestamp
			//cout << "The content of the string data is: " << bufferobj << endl;
			strcpy(datastr,bufferobj.c_str());
			//cout << "The content of the string data is: " << datastr << endl;
			
			strncpy(month,datastr,2);
			strncpy(day,(datastr+3),2);
			strncpy(year,(datastr+6),4);
			strncpy(hour,(datastr+11),2);
			strncpy(min,(datastr+14),2);
			strncpy(sec,(datastr+17),2);
			td.Set(atoi(year),atoi(month),atoi(day),atoi(hour),atoi(min),atoi(sec));
			//cout << month << endl;
			//cout << day << endl;
			//cout << year << endl;
			//cout << hour << endl;
			//cout << min << endl;
			//cout << sec << endl;
			unixtime = td.Convert(); //Here the unix timestamp is set to the acquisition starting time in LTC
			//cout << unixtime << endl;
			
		}
		
		if (i==9){//Acquisition live time (the first of the two numbers)
			bufferobj = bufferobj.substr(0,bufferobj.find(" "));
			//cout << "Acquisition livetime line: " << buffer << endl;
			//cout << "Acquisition livetime line: " << bufferobj << endl;
			//spe >> buffer;
			//cout << "Acquisition livetime string: " << buffer << endl;
			time = atof(bufferobj.c_str());
			//time = (Double_t)atof(bufferstr.c_str());
			cout << "Acquisition livetime: " << time << " secs" << endl;
		}
		
		if (i==11){//Info on how many channels are acquired (the second of the two numbers)
			//cout << "Channels line: " << bufferobj.substr(bufferobj.find(" ")+1,bufferobj.find("\n")) << endl;
			bufferobj = bufferobj.substr(bufferobj.find(" ")+1,bufferobj.find("\n"));
			//spe >> channelstr;//The first value is the number of first channel
			//spe >> channelstr;//The second value is the number of last channel
			channels = atoi(bufferobj.c_str());
			channels++;//The number of channels is 1 more of the last channel number
			cout << "Number of channels: " << channels << endl;
			
			break;
		}
		
		
		//spe.getline(buffer,sizeof(buffer));
	}
	
	
	//Loading the ADC counts from the file and put them in the vector entry_tmp
	while (spe.good()) {
		//The only content of each line is the number of counts in the channel
		getline(spe,bufferobj);
		entry = atof(bufferobj.c_str());
		entries -> push_back(entry);
		if (entries -> size() >= channels) break;
		if (!spe.good()) break;
	}
	
	spe.close();
		
	
	//cout << "Number of vector slots: " << entries->size() << endl;
	
	//Create and fill the histogram
	TH1D* histoADC = new TH1D("histoADC","",entries->size(),1,entries->size());
	for(Int_t i=0; i<entries->size(); i++){
		histoADC->Fill(i+1,entries->at(i));
	}
	
	delete entries;
	
	return histoADC;

}
//------------------------------------------//

#if defined(__CLING__)
TH1D* Gator::loadSpe(const char* dir, Double_t& aqtime)
#else
TH1D* loadSpe(const char* dir, Double_t& aqtime)
#endif
{

	//-------------------------------------------------//
	// Load of the sample histogram from the SPE files //
	// this version is only for calibration files      //
	//-------------------------------------------------//

	string tmpstr("");
	stringstream tmpsstr;
	
	aqtime=0;
	
	Double_t time_tmp=0, entry=0; 
	Int_t channel=0;
	
	vector<Double_t>* entries_tmp = new vector<Double_t>;
	vector<Double_t>* entries = new vector<Double_t>;
	
	string dirname(dir);
	string command("ls ");
	vector<string> spelist;
	
	command += dirname;
	command += string("*.Spe > tmplist");
	
	system(command.c_str());
	ifstream tmplist("tmplist");
	
	while(getline(tmplist,tmpstr)){
		spelist.push_back(tmpstr);
	}
	
	int nSpe = spelist.size();
	
	system("rm -f tmplist");
	 
	ifstream spe;
	
	
	
	for(int iSpe=0; iSpe<nSpe; iSpe++){
		
		string infilename = spelist.at(iSpe);
		cout << "\nLoading file <" << infilename << ">" << endl;
		
		spe.open(infilename.c_str());
		
		if (spe.fail()){
			cout << "Error: can't open the file <" << infilename.c_str() << ">\nSkipping the entry!	" << endl;
			continue;//Skip the loop
		}
		
		entries_tmp->clear();
		
		int channels;
		
		for (Int_t line = 0; line<=11; line++) {
			//Here I go ahead in the file discarding the header lines
			getline(spe,tmpstr);
			
			if (line==9){//Acquisition live time (the first of the two numbers)
				//cout << "Content of \"tmpstr\" at line 9 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				time_tmp = atof(tmpstr.c_str());
				cout << "Acquisition time: " << time_tmp << " secs" << endl;
			}
			
			if (line==11){//Info on how many channels are acquired (the second of the two numbers)
				//cout << "Content of \"tmpstr\" at line 11 is: " << tmpstr << endl;
				tmpsstr.str("");
				tmpsstr << tmpstr;
				tmpsstr >> tmpstr;
				tmpsstr >> tmpstr;
				channels = atoi(tmpstr.c_str());
				channels++;//The number of channels is 1 more of the last channel number
				cout << "Number of channels: " << channels << endl;
			}
		}
		
		//Loading the ADC counts from the file and put them in the vector entry_tmp
		while (spe.good()) {
			spe >> tmpstr;
			entry = atof(tmpstr.c_str());
			if(iSpe==0){
				entries -> push_back(entry);
				if (entries->size()==channels) break;
			} else {
				entries_tmp -> push_back(entry);
				if (entries_tmp->size()==channels) break;
			}
				
			if (!spe.good()) break;
		}
		
		if(iSpe==0){
			cout << "Data vector size: " << entries->size() << endl;
		}else{
			cout << "Data vector size: " << entries_tmp->size() << endl;
			for(Int_t i=0; i<entries_tmp->size(); i++){
				(*entries)[i] += entries_tmp->at(i);
			}
		}
		
		spe.close();
		
		aqtime += time_tmp; //Stored total total acquisition time
		
	}
	
	//Create and fill the histogram
	TH1D* histoMCA = new TH1D("histoMCA","Spectrum ADC",entries->size(),1,entries_tmp->size());
	for(Int_t i=0; i<entries->size(); i++){
		histoMCA->Fill(i+1,entries->at(i));
	}
	
	//-------------------------------------------------------------//	
	// Finished to load of the sample histogram from the SPE files //
	//-------------------------------------------------------------//

	
	delete entries_tmp;
	delete entries;
	
	return histoMCA;

}
//------------------------------------------//

#if defined(__CLING__)
TH1D* Gator::loadSingleSPE(const char* name, Double_t& aqtime)
#else
TH1D* loadSingleSPE(const char* name, Double_t& aqtime)
#endif
{
	
	Int_t nchans;
	ifstream spe;
	string tmp_str;
	cout << "\nOpening file <" << name << ">" << endl;
	
	spe.open(name);
	
	if(spe.fail()){
		cout << "Error: can't open the file <" << name << ">\nReturning 0 as TH1F pointer!" << endl;
		return 0;//Skip the loop
	}
	
	stringstream ss_tmp;
	
	//Getting rid of the first 11 lines (the header)
	for(int i = 0; i<=11; i++){
	
		getline(spe,tmp_str);
		
		if (i==9){//Acquisition lifetime
			ss_tmp.str("");
			ss_tmp << tmp_str;
			ss_tmp >> tmp_str;
			aqtime = atof(tmp_str.c_str());
		}
		
		if (i==11){//Taking the number of MCA channels
			ss_tmp.str("");
			ss_tmp << tmp_str;
			ss_tmp >> tmp_str;
			ss_tmp >> tmp_str;//The second column is the last channel
			nchans = atoi(tmp_str.c_str()) + 1;
		}
	}
	
	
	TH1D* histoMCA = new TH1D("histoMCA","Spectrum MCA",nchans,1,nchans);
	
	//Start getting the channels entries
	for(int i=0; i<nchans; i++){
		spe >> tmp_str;
		histoMCA->Fill(i+1,atof(tmp_str.c_str()));
	}
	
	spe.close();
	
	return histoMCA;
	
}


