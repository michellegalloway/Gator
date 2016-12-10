#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;


//This is for a sharp transfer of sample from exposure to "underground"
double SatActivityA(double activ_days, double cool_days, double half_life, double A0){
	
	const Double_t day = 24.*3600.; //Secs in a day
	//t_cool = cool_days*day; //cool down time (in secs)
	//Double_t activ_time = 348.; //In days
	
	//Double_t A0 = 0.12; //At the start of the measure (in mBq/kg)
	//Double_t half_life = 1925.28; //In days
	
	
	//Convert everything in seconds
	//activ_days *= day;
	//half_life *= day;
	
	Double_t tau = half_life/log(2);
	Double_t A_sat = A0*exp(cool_days/tau)/( 1-exp(-activ_days/tau) );
	
	cout << "Saturation activity: " << A_sat << " mBq/kg" << endl << endl;
	
	return A_sat;
	
}

//This is for a trasition from the full activation to a cool down passing through a "nuisance activation"" period. The nuisance activation strenght is given as the ratio to the full activation strenght
double SatActivityB(double activ_days, double cool_days, double nuis_days, double activ_ratio, double half_life, double A0){
	
	double tau = half_life/log(2);
	double den = 1.-exp(-activ_days/tau) + activ_ratio*( exp(nuis_days/tau)-1. );
	double A_sat = A0*exp((cool_days+nuis_days)/tau)/den;
	
	cout << "Saturation activity: " << A_sat << " mBq/kg" << endl << endl;
	
	return A_sat;
	
}