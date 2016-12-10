#include "trigrate.hh"
//#include "SPEManager.hh"

using namespace std;

trigrate::trigrate():ch_date1(TDatime(2011,8,16,0,0,0))
{
	
	m_filename = "";
	
	m_workdir = "./";
	
	m_ptimestamp  = new std::vector<ULong_t>();
	m_ptrigrate = new std::vector<Double_t>();
	m_ptrigrate_err = new std::vector<Double_t>();
	
	//m_pSPEmanager = 0;
	
	m_sorted = false;
	m_txtfileflg = false; //No data loaded.... it will be taken from SPE file
	
	m_newtimestamp = new std::vector<ULong_t>();
	m_newtrigrate = new std::vector<Double_t>();
	m_newtrigrate_err = new std::vector<Double_t>();
	
}


trigrate::trigrate(const Char_t* filename, const Char_t* workdir):ch_date1(TDatime(2011,8,16,0,0,0))
{
	
	m_workdir = workdir;
	if(m_workdir.substr(m_workdir.size()-1,1)!=string("/")) m_workdir = m_workdir + "/";
	m_filename = filename;
	m_ptimestamp  = new std::vector<ULong_t>();
	m_ptrigrate = new std::vector<Double_t>();
	m_ptrigrate_err = new std::vector<Double_t>();
	
	//m_pSPEmanager = 0;
	
	m_txtfileflg = false; //Loading data from txt file like "triggerrate.txt"
	m_sorted = false;
	
	//cout << "trigrate::constructor(): Calling the routine LoadTxt()." << endl;
	LoadTxt();
	
	m_newtimestamp = new std::vector<ULong_t>();
	m_newtrigrate = new std::vector<Double_t>();
	m_newtrigrate_err = new std::vector<Double_t>();

}


trigrate::~trigrate()
{
	
	if(m_ptimestamp) delete m_ptimestamp;
	if(m_ptrigrate) delete m_ptrigrate;
	if(m_ptrigrate_err) delete m_ptrigrate_err;
	if(m_newtimestamp) delete m_newtimestamp;
	if(m_newtrigrate) delete m_newtrigrate;
	if(m_newtrigrate_err) delete m_newtrigrate_err;
	//if(m_pSPEmanager) delete m_pSPEmanager;
}


Double_t trigrate::GetEntry(Int_t i)
{
	
	if(!m_txtfileflg){
		cout<<"\ntrigrate::GetEntry(...): No data from static archive was loaded returning -1!!!\n\n"<<endl;
		return -1;
	}
	
	if(!m_sorted) SortData();
	
	return m_ptrigrate->at(i);
	
}


Double_t trigrate::GetEntryErr(Int_t i)
{
	
	if(!m_txtfileflg){
		cout<<"\ntrigrate::GetEntryErr(...): No data from static archive was loaded returning -1!!!\n\n"<<endl;
		return -1;
	}
	
	if(!m_sorted) SortData();
	
	return m_ptrigrate_err->at(i);
	
}


void trigrate::SetWorkDir(const Char_t* dir)
{
	
	m_workdir = dir;
	if(m_workdir.substr(m_workdir.size()-1,1)!=string("/")) m_workdir = m_workdir + "/";
	
	return;
}


void trigrate::LoadTxt()
{
	cout << "trigrate::LoadTxt(): Start of routine." << endl;
	ifstream infile(m_filename.c_str());
	
	std::string str_timestamp, str_trigrate, str_trigrate_err;
	ULong_t timestamp;
	Double_t trigrate,trigrate_err;
	
	//cout << "trigrate::LoadTxt(): Entering 'while loop'." << endl;
	while(infile.good()){
		
		infile >> str_timestamp >> str_trigrate >> str_trigrate_err;
		
		if(!infile.good()) break;
		
		timestamp = (ULong_t)atof(str_timestamp.c_str());
		trigrate = atof(str_trigrate.c_str());
		trigrate_err = atof(str_trigrate_err.c_str());
		
		m_ptimestamp -> push_back(timestamp);
		m_ptrigrate -> push_back(trigrate);
		m_ptrigrate_err -> push_back(trigrate_err);
		
	}
	
	infile.close();
	
	m_datasize = m_ptimestamp->size();
	
	for(int i=0; i<m_datasize; i++){
		
		if(i==0){
			m_minTstamp=m_ptimestamp->at(i);
			m_maxTstamp=m_ptimestamp->at(i);
		}else{
			if(m_ptimestamp->at(i)<m_minTstamp) m_minTstamp = m_ptimestamp->at(i);
			if(m_ptimestamp->at(i)>m_maxTstamp) m_maxTstamp = m_ptimestamp->at(i);
		}
		
	}
	
	cout << "trigrate::LoadTxt(): End of routine. Data points: " << m_datasize << endl;
	
	m_txtfileflg = true;
	
	return;	
}


void FindSPE(const Char_t* dirname=0, Bool_t update_flg=false, Bool_t onlysource_flg=false)
{
	
	

}


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
			
			s_command = string("ls ") + s_dirname + string("sources/*.Spe >> tmplist");
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


void trigrate::FindSPE(const Char_t* dirname, Bool_t update_flg, Bool_t onlysource_flg)
{
	
	//This is the function how was defined in the very first version of the "trigrate" class
	
	FindSPE(dirname, (double)545, (double)16380, update_flg, onlysource_flg);
	
	return;
	
}


TH1F* trigrate::LoadSPE(const Char_t* filename,
				Double_t& time,
				ULong_t& unixtime,
				string* descr_str)
{//This pointer can be set to 0 if description is not wanted
	

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



void trigrate::SortData(){
	
	cout << "trigrate::SortData(): Start of routine." << endl;
	
	ULong_t tmp_timestamp;
	Double_t tmp_trigrate;
	Double_t tmp_trigrate_err;
	
	Bool_t sorted = false;
	
	
	for (Int_t i=0; i<m_datasize; i++) {
		if(sorted==true) return;
		sorted = true;
		for (Int_t j=m_datasize-1; j > i; j--){
			if((*m_ptimestamp)[j] < (*m_ptimestamp)[j-1]){
				sorted = false;
				
				tmp_timestamp = (*m_ptimestamp)[j-1];
				(*m_ptimestamp)[j-1] = (*m_ptimestamp)[j];
				(*m_ptimestamp)[j] = tmp_timestamp;
				
				tmp_trigrate = (*m_ptrigrate)[j-1];
				(*m_ptrigrate)[j-1] = (*m_ptrigrate)[j];
				(*m_ptrigrate)[j] = tmp_trigrate;
				
				tmp_trigrate_err = (*m_ptrigrate_err)[j-1];
				(*m_ptrigrate_err)[j-1] = (*m_ptrigrate_err)[j];
				(*m_ptrigrate_err)[j] = tmp_trigrate_err;
				
			}
		}
	}
	
	
	m_sorted = true;
	cout << "trigrate::SortData(): End of routine." << endl;
	return;
}



TGraphErrors* trigrate::MakeFullPlot()
{
	
	cout << "trigrate::MakeFullPlot(...): Start of routine." << endl;
	
	if((m_ptimestamp==0) || (m_ptrigrate==0) || (m_ptrigrate_err==0)){
		cout << "trigrate::MakeFullPlot(): data is not loaded. At least one of the 3 data vectors is not allocated.\nNo plot is going to be generated!\n\n" << endl;
		return 0;
	}
	

	TGraphErrors* pgraph = new TGraphErrors();
	
	Int_t point = 0;
	
	for(int i=0; i<m_ptimestamp->size(); i++){
		
		pgraph -> SetPoint(point,(Double_t)m_ptimestamp->at(i),m_ptrigrate->at(i));
		pgraph -> SetPointError(point,0,m_ptrigrate_err->at(i));
		point++;
	}
	
	cout << "trigrate::MakeFullPlot(...): End of routine." << endl;

	return pgraph;
}



TGraphErrors* trigrate::MakeGraph(ULong_t startime, ULong_t endtime)
{

	cout << "trigrate::MakeGraph(...): Start of routine." << endl;

	if((m_ptimestamp==0) || (m_ptrigrate==0) || (m_ptrigrate_err==0)){
		cout << "trigrate::MakeGraph(...): data is not loaded. At least one of the 3 data vectors is not allocated.\nNo plot is going to be generated!\n\n" << endl;
		return 0;
	}
	

	TGraphErrors* pgraph = new TGraphErrors();
	
	Int_t point = 0;
	for(int i=0; i<m_ptimestamp->size(); i++){
		if((startime<=m_ptimestamp->at(i)) && (m_ptimestamp->at(i)<=endtime)){
		
			pgraph -> SetPoint(point,(Double_t)m_ptimestamp->at(i),m_ptrigrate->at(i));
			pgraph -> SetPointError(point,0,m_ptrigrate_err->at(i));
			point++;
		}
	}
	
	cout << "trigrate::MakeGraph(...): End of routine." << endl;

	return pgraph;
}



void trigrate::WriteData(const Char_t* filename, Bool_t append)
{
	
	ofstream newtrigfile;
	
	if(append){//Append to the file only the new found timestamps
		newtrigfile.open(filename,ios_base::app);
		for(int i=0; i<m_newtimestamp->size(); i++){
			newtrigfile << m_newtimestamp->at(i) << "\t" << m_newtrigrate->at(i) << "\t" << m_newtrigrate_err->at(i) << endl;
		}
		
	}else{//Overwriting the file
		newtrigfile.open(filename);
		for(Int_t i=0; i<m_datasize; i++){
			newtrigfile << m_ptimestamp->at(i) << "\t" << m_ptrigrate->at(i) << "\t" << m_ptrigrate_err->at(i) << endl;
		}
	}
	
	
	
	newtrigfile.close();
	
}



void trigrate::PrintDates(Int_t num)
{
	
	if(!m_sorted) SortData();
	if(num > m_datasize){
		num = m_datasize; //protection!!!
	}
	
	cout << "The last " << num << " trigger aquisition of " << m_datasize << " are:" << endl;
	
	TDatime t;
	for(int i=m_datasize-num; i<m_datasize; i++){
		t.Set((*m_ptimestamp)[i]);
		cout << (*m_ptimestamp)[i] << "  -->  " << t.GetYear() << " " << t.GetMonth() << " " << t.GetDay() << "   " << t.GetHour() << ":" << t.GetMinute() << endl;
	}
}



void trigrate::RepairTxt(const Char_t* SPEdir)
{
	
	if(!m_txtfileflg){
		cout << "\ntrigrate::RepairTxt(...): There is not any loaded static archive loaded for the trigger rate. Termitatio of the routine withot doing anything.\n\n" << endl;
	}
	
	
	//Saving the old static archive in a "backup" file
	string backup_comm = string("cp ") + m_filename + string(" ") + m_workdir + m_filename + string(".broken");
	
	system(backup_comm.c_str());
	
	/*
	backup_comm += m_filename;
	backup_comm += " ";
	backup_comm += m_filename.c_str();
	backup_comm += ".broken";
	*/
	
	system(backup_comm.c_str());
	
	SortData();
	FindSPE(SPEdir,false,false); //Called in no update mode (overwriting) and in no "find source only".
	
	
	//Find the min and max timestamp for the data from the new SPE files found by the "FindSPE(...)" routine
	ULong_t mintimestp, maxtimestp;
	
	for(Int_t i=0; i<m_newtimestamp->size(); i++){
		if(i==0){
			mintimestp = m_newtimestamp->at(i);
			maxtimestp = mintimestp;
		}else{
			if(mintimestp > m_newtimestamp->at(i)) mintimestp = m_newtimestamp->at(i);
			if(maxtimestp < m_newtimestamp->at(i)) maxtimestp = m_newtimestamp->at(i);
		}
	}
	
	
	//Save in the tmp vectors only the data that have timestamp < mintimestamp
	std::vector<ULong_t>* tmp_timestamp = new std::vector<ULong_t>();
	std::vector<Double_t>* tmp_trigrate = new std::vector<Double_t>();
	std::vector<Double_t>* tmp_trigrate_err = new std::vector<Double_t>();
	
	for(Int_t i=0; i<m_datasize; i++){
		if((m_ptimestamp->at(i)<mintimestp)){
			tmp_timestamp->push_back(m_ptimestamp->at(i));
			tmp_trigrate->push_back(m_ptrigrate->at(i));
			tmp_trigrate_err->push_back(m_ptrigrate_err->at(i));
		}
	}
	
	
	//The normal vectors will have only the data old data before the mintimestamp
	(*m_ptimestamp) = (*tmp_timestamp);
	(*m_ptrigrate) = (*tmp_trigrate);
	(*m_ptrigrate_err) = (*tmp_trigrate_err);
	
	
	//To the data vectors will be add new data found by the routine
	for(Int_t i=0; i<m_newtimestamp->size(); i++){
		m_ptimestamp->push_back(m_newtimestamp->at(i));
		m_ptrigrate->push_back(m_newtrigrate->at(i));
		m_ptrigrate_err->push_back(m_newtrigrate_err->at(i));
	}
	
	m_datasize = m_ptimestamp->size();
	
	SortData();
	
	WriteData(m_filename.c_str(), false); //Overwriting the data.
	

	if(tmp_timestamp) delete tmp_timestamp;
	if(tmp_trigrate) delete tmp_trigrate;
	if(tmp_trigrate_err) delete tmp_trigrate_err;

	
	return;
}


/*
void trigrate::UpdateTrigRate(Bool_t source)
{
	
	if(m_pSPEmanager==0){
		cout << "\ntrigrate::UpdateTrigRate(...): SPEManager not set inside the class.\nThis method is thought to be used only by an SPEManager class!!!\n\n" << endl;
		return;
	}
	
	m_filename = m_pSPEmanager->GetStaticArchive(); //I have now the name with the path of the "triggerrate.txt" static archive
	
	string workdir = m_pSPEmanager->GetWorkDir();
	
	vector<string>* newspefiles;
	if(!source) {
		newspefiles = m_pSPEmanager->GetNewSPEVector();
	}else{
		newspefiles = m_pSPEmanager->GetNewSourceVector();
	}
	
	TH1F* histo = 0;
	
	m_newtimestamp->clear();
	m_newtrigrate->clear();
	m_newtrigrate_err->clear();
	
	Double_t aqtime;
	ULong_t startimstamp;
	
	Int_t newdatasize = newspefiles->size();
	
	for(int i=0; i<newdatasize; i++){
		
		Double_t rate, rate_err;
		
		string spename = workdir+newspefiles->at(i);
		
		if(histo) delete histo;
		
		histo = LoadSPE(spename.c_str(),aqtime,startimstamp);
		
		if(!histo){
			cout << "\n\ntrigrate::UpdateTrigRate(): File <"<< spename <<"> have not given a valid histogram. Update of the static archive will not happen!!!\nTry to use manually the spefind or the repair routines.\n\n" << endl;
			return;
		}
		
		
		rate = (histo->Integral(545,16380))/aqtime; //Hard-coded be careful!!!
		rate_err = TMath::Sqrt(histo->Integral(545,16380))/aqtime; //Hard-coded be careful!!!
		
		//cout << startimstamp << "\t" << rate << "\t" << rate_err << endl;
		//if(!listfile.good()) break;
		
		cout << "Start timestamp: " << startimstamp << endl;
		cout << "DAQ rate: (" << rate << " +- " << rate_err << ") secs^{-1}" << endl;
		
		m_newtimestamp->push_back(startimstamp);
		m_newtrigrate->push_back(rate);
		m_newtrigrate_err->push_back(rate_err);
		
	}
	
	if((newdatasize==m_newtimestamp->size())&&(newdatasize==m_newtrigrate->size())&&(newdatasize==m_newtimestamp->size())){
		cout << "\n\ntrigrate::UpdateTrigRate(): There are " << newdatasize << " new entries to add at the static archive.\n\n" << endl;
		WriteData(m_filename.c_str(),true);
	} else {
		cout << "\n\ntrigrate::UpdateTrigRate(): Data cannot ba updated because of mismatching between vectors dimension!!!\n\n" << endl;
	}
	
	
	return;
}
*/