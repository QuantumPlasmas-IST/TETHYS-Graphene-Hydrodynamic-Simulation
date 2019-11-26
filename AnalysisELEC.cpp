#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <ctime>
#include <numeric>
#include <fftw3.h>

#include "dyakonovshur.h"


using namespace std;




int main(int argc, char **argv){

	cout << "\n";
	cout << "*******************************************************"<<endl;
	cout << "****************** Antenna Analysis *******************"<<endl;
	cout << "*******************************************************"<<endl;	
	

	int N = atoi(argv[1]);	
	
	float S = atof(argv[3]);	
	
	float *Time;			
	Time =(float*) calloc (N,sizeof(float));
	
	float *in_potential;			
	in_potential =(float*) calloc (N,sizeof(float));
	
	float *in_current;			
	in_current =(float*) calloc (N,sizeof(float));
	
	float *in_dipole;			
	in_dipole =(float*) calloc (N,sizeof(float));
	
	float *in_D_dipole;			
	in_D_dipole =(float*) calloc (N,sizeof(float));


	float *in_K_energy;			
	in_K_energy =(float*) calloc (N,sizeof(float));


	float *out_D2_dipole;			
	out_D2_dipole =(float*) calloc (N,sizeof(float));

	float *out_D2_dipoleAux;			
	out_D2_dipoleAux =(float*) calloc (N,sizeof(float));


	float *out_power;			
	out_power =(float*) calloc (N,sizeof(float));

	//////////////////Reading input file ///////////////////////////////
	ifstream input;
	input.open(argv[2]);

	cout << "Reading input file";
	int i=0;
	if(input.is_open())
	{
		while(input.good())
		{
			input >> Time[i] >> in_potential[i] >> in_current[i] >>in_dipole[i] >> in_D_dipole[i] >> in_K_energy[i];
			i++;	
		}
	}
	cout << "\tDONE!"<<"\tS="<<S<<endl;
	float dt = Time[2]- Time[1];
	cout << "dt " << dt <<endl;
	cout << "time max " << Time[N-1] <<endl;	
	////////////////////////////////////////////////////////////////////


	ofstream logfile;
	logfile.open("Antenna.log",std::ios_base::app);
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	logfile << "\n#Simulation @ " << asctime(timeinfo) ;
	
	
	///// Gaussian convolution to obtain dipole acceleration ///////////
	float t = 65;
	int M = 50;
    ConvolveGauss(1, M, t, in_D_dipole, out_D2_dipole, N);
	////////////////////////////////////////////////////////////////////
	
	
	/////  obtain dipole acceleration ///////////
	for(int k=1;k<N-1;k++){
	
	out_D2_dipoleAux[k] = (-0.5*in_D_dipole[k-1]+0.5*in_D_dipole[k+1])/dt;
	
	}
    
	////////////////////////////////////////////////////////////////////
	

	//////////////// Writing output file ///////////////////////////////
	string antenafile = "antena_" + to_string(S)+ ".dat" ;
	antenafile.erase (antenafile.end()-9, antenafile.end()-5);
	ofstream data_Antena;
	data_Antena.open(antenafile);
	data_Antena << scientific; 
	
	for(int k=0;k<N;k++){
		out_power[k] = out_D2_dipole[k]*out_D2_dipole[k];
		data_Antena <<Time[k] <<"\t"<<in_potential[k] <<"\t"<< in_current[k] <<"\t"<< in_K_energy[k] <<"\t"<< out_D2_dipole[k]/Time[N-1] <<"\t"<< out_power[k]/(Time[N-1]*Time[N-1])  <<endl; 
	}
	////////////////////////////////////////////////////////////////////



	///// RMS values  //////////////////////////////////////////////////
	
	float I_rms, Ug_rms, d2P_rms,Power_rms; 
	Ug_rms = RMS(N, dt, in_potential);
	I_rms = RMS(N, dt, in_current);
	d2P_rms = RMS(N, dt, out_D2_dipole);
	Power_rms = RMS(N, dt, out_power);
	
	cout << "RMS Ugs : " << Ug_rms <<endl;
	cout << "RMS Ids : " << I_rms  <<endl;
	cout << "RMS accel dipole : " << d2P_rms  <<endl;
	cout << "RMS power : " << Power_rms  <<endl;
	////////////////////////////////////////////////////////////////////
	logfile << "#S \t V_gs \t I_ds \t d2P/dt2 \t Power_rms" <<endl;
	logfile << S <<"\t"<< Ug_rms <<"\t"<< I_rms <<"\t"<< d2P_rms <<"\t"<< Power_rms<<endl;

	data_Antena.close();
	input.close();



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Power Analysis 
//
////////////////////////////////////////////////////////////////////////	
	
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Current Analysis
//
////////////////////////////////////////////////////////////////////////	
	
	
	//	 Extrema and Growth Rate
	string extremafileI = "current_extrema_" + to_string(S)+ ".dat" ;
	extremafileI.erase (extremafileI.end()-9, extremafileI.end()-5);
	float saturation = 0.0;
	float tau = 0.0;
	float error = 0.0;
	ExtremaFinding(in_current, N, S, dt,saturation,tau, error, extremafileI);

	cout << "Current saturation value: " << saturation << endl;
	cout << "Time for 99% of saturation: " << tau <<endl;			

	logfile << "#S \t I_sat \t Time for 99% of saturation \t time error" <<endl;
	logfile << S  <<"\t"<< saturation <<"\t"<< tau <<"\t"<< error <<endl;
	
	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// VDS Analysis
//
////////////////////////////////////////////////////////////////////////	
	
	
	//	 Extrema and Growth Rate
	string extremafileV = "tension_extrema_" + to_string(S)+ ".dat" ;
	extremafileV.erase (extremafileV.end()-9, extremafileV.end()-5);
	saturation = 0.0;
	tau = 0.0;
	error = 0.0;
	ExtremaFinding(in_potential, N, S, dt,saturation,tau, error, extremafileV);

	cout << "Tension saturation value: " << saturation << endl;
	cout << "Time for 99% of saturation: " << tau <<endl;			

	logfile << "#S \t VDS_sat \t Time for 99% of saturation \t time error" <<endl;
	logfile << S  <<"\t"<< saturation <<"\t"<< tau <<"\t"<< error <<endl;	
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Redimension 
//
////////////////////////////////////////////////////////////////////////	

	string satantenafile = "sat_antena_" + to_string(S)+ ".dat" ;
	satantenafile.erase (satantenafile.end()-9, satantenafile.end()-5);
	ofstream data_Saturation;
	data_Saturation.open(satantenafile);
	data_Saturation << scientific; 
	
	int cut;
	cut = tau/dt;	
	
	
	
	for(int k=0;k<N;k++){
		if(k>cut){
			data_Saturation <<Time[k] <<"\t"<< in_potential[k]<<"\t"<< in_current[k] <<"\t"<< in_K_energy[k] <<"\t"<< out_D2_dipole[k] <<"\t" << out_power[k] <<endl; 
		}
	}

	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	logfile.close();
	cout << "*******************************************************"<<endl;	
	return 0;
	
}
