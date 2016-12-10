#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

void ScaleMeanActiv(double meanactiv,double runtime,double halflife){
	
	/*
	double meanact = 0.76;
	double runtime = 4.;
	double halflife = 1925.28;
	*/
	
	double tau = halflife/log(2);
	double init_act = meanactiv*runtime/( tau*( 1.-exp(-runtime/tau)) );
	
	cout << "Rescaled activity: " << init_act << " mBq/kg" << endl << endl;
	
	return;

}