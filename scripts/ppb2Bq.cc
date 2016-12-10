#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
using namespace std;


namespace isotopes{
}


//OBSOLETE
double ppb2Bq(const char element[] = "U238", double ppb = 1e-3) 
{
  string material = element;
  
  double A, T12;
  
  if(material == "U238")
    { A = 238; T12 = 4.468e+9 * 365.*24.*3600.; }
  else if(material == "Ra226")
    { A = 226; T12 = 1600.0 * 365.*24.*3600.; }
  else if(material == "U235")
    { A = 235; T12 = 7.038e+8 * 365.* 24.*3600.; }
  else if(material == "Th232")
    { A = 232; T12 = 1.405e10 * 365.*24.*3600.; }
  else if(material == "Th228")
    { A = 228; T12 = 1.9131 * 365.*24.*3600.; }
  else if(material == "K40")
    { A = 40; T12 = 1.248e+9 * 365.*24.*3600.; }
  else if(material == "Pb210")
    { A = 210; T12 = 22.3 * 365.*24.*3600.; }
  else if(material == "Kr85")
    { A = 85; T12 = 10.756 * 365.*24.*3600.; }
  else if(material == "Co60")
    { A = 60; T12 = 1925.28*24.*3600.; }//Lifetime directly in days (from Nudat)
  else if(material == "Cs137")
    { A = 137; T12 = 30.08*365.*24.*3600.; }//Lifetime directly in days (from Nudat)
  else{
	  return -1;
  }
  
  // lifetime in seconds
  double tau = T12/log(2);
  // Avogadro number
  double Na = 6.022e23; // mol-1
  
  double conc  = 1e-9; // ppb

  // COMPUTE
  double R_g 	 = (conc*Na)/(A*tau); 
  double R_kg    = R_g*1000; 
  
  cout << endl;
  cout << material << " --> Activity of " << conc << " = " << R_kg << " Bq/kg " << endl;
  double ScalingFactor = 1/R_kg;
  cout << "  1 Bq/kg = " << ScalingFactor << " ppb " << endl;
  cout << endl;
  cout << " ==> For " << ppb << " ppb => " << ppb*R_kg*1000. << " mBq/kg " << endl;
  cout << endl;
	
  return ppb*R_kg*1000.;
	
}

double ppb2mBq(const char element[] = "U238", double ppb = 1.) 
{
  string material = element;
  
  double A, T12;
  
  if(material == "U238")
    { A = 238; T12 = 4.468e+9 * 365.*24.*3600.; }
  else if(material == "Ra226")
    { A = 226; T12 = 1600.0 * 365.*24.*3600.; }
  else if(material == "U235")
    { A = 235; T12 = 7.038e+8 * 365.* 24.*3600.; }
  else if(material == "Th232")
    { A = 232; T12 = 1.405e10 * 365.*24.*3600.; }
  else if(material == "Th228")
    { A = 228; T12 = 1.9131 * 365.*24.*3600.; }
  else if(material == "K40")
    { A = 40; T12 = 1.248e+9 * 365.*24.*3600.; }
  else if(material == "Pb210")
    { A = 210; T12 = 22.3 * 365.*24.*3600.; }
  else if(material == "Kr85")
    { A = 85; T12 = 10.756 * 365.*24.*3600.; }
  else if(material == "Co60")
    { A = 60; T12 = 1925.28*24.*3600.; }//Lifetime directly in days (from Nudat)
  else if(material == "Cs137")
    { A = 137; T12 = 30.08*365.*24.*3600.; }//Lifetime directly in days (from Nudat)
  else if(material == "Mn54")
    { A = 54; T12 = 312.12*365.*24.*3600.; }//Lifetime directly in days (from Nudat)
  else{
	  return -1;
  }
  
  // lifetime in seconds
  double tau = T12/log(2);
  // Avogadro number
  double Na = 6.022e23; // mol-1
  
  double Ratio = (1e-3*Na)/(A*tau); //Ratio between activity in mBq/kg and concentration in ppb

  // COMPUTE
  double Activity = ppb*Ratio;
  
  cout << endl;
  cout << material << ": --> Activity of 1 ppb = " << Ratio << " mBq/kg " << endl;
  cout << " ==> For " << ppb << " ppb " <<  "--> " << Activity << " mBq/kg " << endl << endl;
	
  return Activity;
	
}


double mBq2ppb(const char element[] = "U238", double Activity=1.){
	
    string material = element;
  
    double A, T12;
  
    if(material == "U238")
      { A = 238; T12 = 4.468e+9 * 365.*24.*3600.; }
    else if(material == "Ra226")
      { A = 226; T12 = 1600.0 * 365.*24.*3600.; }
    else if(material == "U235")
      { A = 235; T12 = 7.038e+8 * 365.* 24.*3600.; }
    else if(material == "Th232")
      { A = 232; T12 = 1.405e10 * 365.*24.*3600.; }
    else if(material == "Th228")
      { A = 228; T12 = 1.9131 * 365.*24.*3600.; }
    else if(material == "K40")
      { A = 40; T12 = 1.248e+9 * 365.*24.*3600.; }
    else if(material == "Pb210")
      { A = 210; T12 = 22.3 * 365.*24.*3600.; }
    else if(material == "Kr85")
      { A = 85; T12 = 10.756 * 365.*24.*3600.; }
    else if(material == "Co60")
      { A = 60; T12 = 1925.28*24.*3600.; }//Lifetime directly in days (from Nudat)
    else if(material == "Cs137")
      { A = 137; T12 = 30.08*365.*24.*3600.; }//Lifetime directly in days (from Nudat)
    else if(material == "Mn54")
      { A = 54; T12 = 312.12*365.*24.*3600.; }//Lifetime directly in days (from Nudat)
	else{
  	  return -1;
    }
	
	
    // lifetime in seconds
    double tau = T12/log(2);
    // Avogadro number
    double Na = 6.022e23; // mol-1
	
	double Ratio = (1e-3*Na)/(A*tau); //Ratio between activity in mBq/kg and concentration in ppb
	
	double ppb = Activity/Ratio;
	
    cout << endl;
    cout << material << ":  --> Concentration of 1 mBq/kg 1 = " << 1./Ratio << " ppb " << endl;
	cout << " ==> For " << Activity << " mBq/kg " <<  "--> " << ppb << " ppb " << endl << endl;
	
	
	return ppb;
	
}
