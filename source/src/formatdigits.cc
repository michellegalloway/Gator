#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;


string formatdigits1(double var, double err, int dig){//By default the digits used for the precision are 2
	
//This function is used to format the output of value with its error with a given number of digits (to be choiced)
	
	ostringstream s_var, s_err;
	ostringstream s_var_sc, s_err_sc;
	string outstring;
	
	int log_v, log_e;
	
	log_e = (int)(log10(err));
	if(log10(err)<0.) log_e--;//In this way works fine always
	
	log_v = (int)(log10(var));
	if(log10(var)<0.) log_v--;//In this way works fine always
	
	int nse = (dig-1) - (int)(log_e);
	if(nse < 0) nse = 0; //I don't need digits at the right of the point
	
	
	s_err_sc.precision(dig-1);
	s_err_sc << scientific << err;
	
	err = atof(s_err_sc.str().c_str());
	
	s_err.precision(nse);
	s_err  << fixed << err;
	//cout << "s_err = " << s_err.str() << endl;
	
	
	//int nse = dig - (int)(log10(err)+1);
	
	//cout << "nse = " << nse << endl;
	
	int nse_sc = log_v - log_e + 1;
	s_var_sc.precision(nse_sc);
	//s_var << setw(nse);
	s_var_sc << scientific << var;
	
	
	var = atof(s_var_sc.str().c_str());
	//err = atof(s_err.str().c_str());
	
	
	//s_err << setprecision(dig) << err;
	s_var.precision(nse);
	s_var << fixed << var;
	//cout << "s_var = " << s_var.str() << endl;
	
	
	
	outstring = s_var.str();
	outstring += " +- ";
	outstring += s_err.str();
	
	return outstring;
}



//------------------------------------------------------------------------------------------------\\

string formatdigits2(double var, int dig){//By default the digits used for the precision are 3
	
//This function is used to format the output of a value with a given number of digits (for example to be used in the upper limits)
	
	ostringstream s_var, s_var_sc;
	
	string outstring;
	
	int log_v;
	
	log_v = (int)(log10(var));
	if(log10(var)<0.) log_v--;//In this way works fine always
	
	int nse = (dig-1) - (int)(log_v);
	if(nse < 0) nse = 0; //I don't need digits at the right of the point
	
	s_var_sc.precision(dig - 1);
	s_var_sc << scientific << var;
	
	var = atof(s_var_sc.str().c_str());
	
	s_var.precision(nse);
	s_var << fixed << var;
	//cout << "s_var = " << s_var.str() << endl;
	
	
	outstring = s_var.str();
	
	return outstring;
}
