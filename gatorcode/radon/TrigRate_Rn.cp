/*
 * TrigRate_Rn.cpp
 *
 */


////////////////////////////////////////////////////////////////////////////////


// Standard libs:
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;


#include "trigrate_v2.hh"
//#include "TrigRate_Rn.hh"
//#include "SPEManager.hh"



////////////////////////////////////////////////////////////////////////////////
// constructor
trigrate::trigrate(): SCVInstMonit()
{
    // blank
}

////////////////////////////////////////////////////////////////////////////////
// destructor
trigrate::~trigrate()
{
    /*    if(m_ptimestamp) delete m_ptimestamp;
     //    if(m_ptrigrate) delete m_ptrigrate;
     //    if(m_ptrigrate_err) delete m_ptrigrate_err;
     if(m_newtimestamp) delete m_newtimestamp;
     if(m_newtrigrate) delete m_newtrigrate;
     if(m_newtrigrate_err) delete m_newtrigrate_err;  */
    //    if(m_pSPEmanager) delete m_pSPEmanager;
}

////////////////////////////////////////////////////////////////////////////////
// define load data function

trigrate::trigrate(string datadir): SCVInstMonit(datadir)
{
    LoadData();
    if(m_dataflag){
		m_size = mv_unixtimes->size();
		SetTmin();
		SetTmax();
        m_YAxisTitle = "Trigger Rate [s^{-1}]";
        m_YAxisRangeMin=0.0001;
        m_YAxisRangeMax=0.01;
	}
}

////////////////////////////////////////////////////////////////////////////////
// load SPE (?) data

void trigrate::LoadData(){
    
    // get files
    string oscommand("ls ");
	oscommand += m_datadir.c_str();
	oscommand += "/*.txt > ./trig.flist";
	system(oscommand.c_str());
    
	ifstream listfiles("./trig.flist");
	string filename;
    
    // read files
    if(!mv_filenames){
        mv_filenames = new vector<string>;
    }else{
        mv_filenames->clear();
    }
    
    if(!mv_values){
        mv_values = new vector<double>;
    }else{
        mv_values->clear();
    }

    if(!mv_errvalues){
        mv_errvalues = new vector<double>;
    }else{
        mv_errvalues->clear();
    }

    if(!mv_unixtimes){
        mv_unixtimes = new vector<time_t>;
    }else{
        mv_unixtimes->clear();
    }
    
    system("rm ./trig.flist");

    
    // read lines
	while (getline(listfiles, filename)) {
        ifstream datafile(filename.c_str());
        string bufferobj; // line with data (first token) and timestamp (second token)
        time_t timestm;  // timestamp
        stringstream tmp_line;
        
        double val;
        double errval;
        int nLine=0;
        
        getline (datafile, bufferobj); //Skip the first line always
        nLine++;
        while (getline (datafile, bufferobj)) {
            nLine++;
            stringstream ss(bufferobj);
            string token;
            int count = 0;
            while(ss>>token){
                count++;
            }
            
            
            // below eliminates lines with not enough values
            if (count != 3) {
                // cout<< "Skipping line " << nLine-1 << " of file <" << filename << ">" << " -> has not enough tokens!" <<endl;
                continue;
            }
            tmp_line.clear();
            tmp_line.str("");
            tmp_line << bufferobj;
            
            tmp_line >> bufferobj; // The first column is the timestamp
            timestm = strtoul(bufferobj.c_str(), NULL, 10); // parses and converts time string to unsigned long integer (base 10)

            
            timestm += 43200;  // //This is to shift the point from the acquisition start time to the acquisition stop time (12 hours). It is much more intuitive for who reads the slow control plots.
            
            if(!tmp_line.good()){
                // cout << "Skipping line " << nLine-1 << " of file <" << filename << ">" << endl;
                continue;
            }
            
            
            tmp_line >> bufferobj; // second column is the value
            val = strtod(bufferobj.c_str(), NULL);
            
            
            // below is hack to eliminate lines with nonsensical (too high value) datapoints
            /*  if (val > 1000){
             cout << "Skipping line " << nLine-1 << " of file <" << filename << ">" << " -> has nonsensical data!" << endl;
             continue;
             } */

            
            tmp_line >> bufferobj; // The third column is the error
            errval = strtod(bufferobj.c_str(), NULL);
            
            // put data filters here if necessary
            
            // get datapoints
            // temp = TempConvert(val);
            mv_unixtimes->push_back(timestm);
            mv_values->push_back(val);
            mv_errvalues->push_back(errval);
            mv_filenames->push_back(filename);
            //            if(listfiles.eof()) break;

        }
        
    }
    
    //Check for consistency of loaded data
    m_dataflag = (mv_values->size()==mv_unixtimes->size())&&(mv_unixtimes->size()==mv_filenames->size());
    // if(m_dataflag) m_size = mv_values->size();

}


/*
////////////////////////////////////////////////////////////////////////////////


void trigrate::FindSPE(const Char_t* dirname, double StartCh, double EndCh, Bool_t update_flg, Bool_t onlysource_flg)
{
	cout << "\ntrigrate::FindSPE: Start of routine" << endl;
	
	m_newtimestamp->clear();
	m_newtrigrate->clear();
	m_newtrigrate_err->clear();
	
	Double_t aqtime, rate, rate_err;
	ULong_t startimstamp;
	TH1F* histo = 0;

	if(StartCh<1) StartCh=1;
	
	string s_dirname = dirname;
	if(s_dirname.substr(s_dirname.size()-1,1)!=string("/")) s_dirname = s_dirname + "/";
	if(!dirname){
		if(!onlysource_flg) system("ls *.SPE > tmplist");
		system("ls *.Spe >> tmplist");
	}else{
		
		string s_command;
		
		if(!onlysource_flg){
			s_command = string("ls ") + s_dirname + string("*.SPE > tmplist");
			//cout << "Command1 = '" << s_command << "'" << endl;
			system(s_command.c_str());
			s_command = string("ls ") + s_dirname + string("sources */  /*   .Spe >> tmplist");

			//cout << "Command2 = '" << s_command << "'" << endl;
			system(s_command.c_str());
		}else{
			s_command = string("ls ") + s_dirname + string("*.Spe >> tmplist");
			//cout << "Command2 = '" << s_command << "'" << endl;
			system(s_command.c_str());
		}
	}
	
	//return;
	
	if((!m_sorted)&&(m_txtfileflg)) SortData();
	
	ifstream listfile("tmplist");
	string spename;
	Bool_t busy_flag = false;
	
	while(listfile.good()){
		busy_flag = false;
		listfile >> spename;
		if(listfile.eof()) break;
		if(histo) delete histo;
		histo = LoadSPE(spename.c_str(),aqtime,startimstamp);
		if(!histo){
			continue;
		}
		
		//Here a starting channel to calculate the spectrum is chosen
		if(ch_date1.Convert()<=startimstamp){
			if(EndCh>16383) EndCh=16383;
			rate = (histo->Integral(StartCh,EndCh))/aqtime;
			rate_err = TMath::Sqrt(histo->Integral(StartCh,EndCh))/aqtime; //Hard-coded be careful!!!
		}else{
			if(EndCh>(12300-11)) EndCh=12300-11;
			rate = (histo->Integral(StartCh,EndCh))/aqtime;
			rate_err = TMath::Sqrt(histo->Integral(StartCh,EndCh))/aqtime; //Hard-coded be careful!!!
		}
		
		//cout << startimstamp << "\t" << rate << "\t" << rate_err << endl;
		//if(!listfile.good()) break;
		
		cout << "Start timestamp: " << startimstamp << endl;
		cout << "DAQ rate: (" << rate << " +- " << rate_err << ") secs" << endl;
		
		//Here I check if this timestamp already exist. 
		//If data is not loaded from a txtfile I skip it.
		if(m_txtfileflg){
			for(Int_t i=0; i<m_datasize; i++){
				if(m_ptimestamp->at(i)==startimstamp){
					busy_flag = true;
					break;
				}
			}
		}
		//cout << m_newtimestamp->size() << "\t" << m_newtrigrate->size() << "\t" << m_newtrigrate_err->size() << endl;
		if(!busy_flag){
			m_newtimestamp->push_back(startimstamp);
			m_newtrigrate->push_back(rate);
			m_newtrigrate_err->push_back(rate_err);
		}
		if(listfile.eof()) break;
	}
	
	
	
	Int_t newdatasize = m_newtimestamp->size();
	cout << "There are " << newdatasize << " new entries.\n" << endl;
	if(histo) delete histo;
	
	//This routine will write only in append mode
	if(update_flg){
		if(newdatasize > 0){
			for(Int_t i=0; i<newdatasize; i++){
				m_ptimestamp->push_back(m_newtimestamp->at(i));
				m_ptrigrate->push_back(m_newtrigrate->at(i));
				m_ptrigrate_err->push_back(m_newtrigrate_err->at(i));
			}
			
			m_datasize = m_ptimestamp->size();
			if(m_txtfileflg){
				WriteData(m_filename.c_str(),true); //Write in append mode
				cout << "Trigger rate file updated with "<< newdatasize << " new entries." << endl;
            }
		}else{
			cout << "No new data in the directory <"<< dirname <<">\n" << endl;
		}
	}
	
	system("rm -f tmplist");
	cout << "trigrate::FindSPE(...): End of routine.\n" << endl;
	return;
}


////////////////////////////////////////////////////////////////////////////////

void trigrate::FindSPE(const Char_t* dirname, Bool_t update_flg, Bool_t onlysource_flg)
{
	//This is the function as it was defined in the very first version of the "trigrate" class
	FindSPE(dirname, (double)545, (double)16380, update_flg, onlysource_flg);
	return;
	
}

 
////////////////////////////////////////////////////////////////////////////////


TH1F* trigrate::LoadSPE(const Char_t* filename,
				Double_t& time,
				ULong_t& unixtime,
				string* descr_str){//This pointer can be set to 0 if description is not wanted
	

	//------------------------------------------------------//	
	// Load of the sample histogram from the SPE files      //
	// This version gives the Unix timestamp and the sample //
	// description if requested                             //
	//------------------------------------------------------//

	std::string bufferobj, data_str("");//, livetimeobj, channelobj;
	
	TDatime td("1970-01-01 00:00:00"); //Used to put in the acquisition data and to generate the correspondent unix timestamp

	//Int_t timeoffset = -3600; //To be used to return the correct Unix timestamp
	
	ifstream spe(filename);
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

			bufferobj = bufferobj.substr(0,bufferobj.find("\n"));
			
			data_str+=bufferobj.substr(6,4); //The year
			data_str+="-";
			data_str+=bufferobj.substr(0,2); //The month
			data_str+="-";
			data_str+=bufferobj.substr(3,2); //The day
			data_str+=" ";
			data_str+=bufferobj.substr(11); //The time (hh:mm:ss)

//			cout << "Acquisition start string: " << data_str << endl;
			td.Set(data_str.c_str());
			cout << "Acquisition start is on: " << td.GetDay() << "/" << td.GetMonth() << "/" << td.GetYear() << "   " << td.GetHour() << ":" << td.GetMinute() << endl;
			td.GetDate();
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
	
	
	Double_t entry;
	std::vector<Double_t>* entries = new std::vector<Double_t>();
	//Loading the ADC counts from the file and put them in the vector entry_tmp
	while (spe.good()) {
		//The only content of each line is the number of counts in the channel
		getline(spe,bufferobj);
		if (!spe.good()) break;
		entry = atof(bufferobj.c_str());
		entries->push_back(entry);
		if (entries->size() >= channels) break;
		if (!spe.good()) break;
	}
	
	spe.close();
		
	
	//cout << "Number of vector slots: " << entries->size() << endl;
	
	//Create and fill the histogram
	TH1F* histoADC = new TH1F("histoADC","",entries->size(),1,entries->size());
	for(Int_t i=0; i<entries->size(); i++){
		histoADC->Fill(i+1,entries->at(i));
	}
	
//	cout << "Deleting vector 'entries'" << endl;
	if(entries) delete entries;
//	cout << "Vector 'entries' deleted" << endl;
	
//	cout << "trigrate::LoadSPE(...) --> End of routine" << endl;	
	return histoADC;
}

 
 */
////////////////////////////////////////////////////////////////////////////////




// The end...
////////////////////////////////////////////////////////////////////////////////

