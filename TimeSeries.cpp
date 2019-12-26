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
//#include <fftw3.h>

#include "dyakonovshur.h"


using namespace std;




int main(int argc, char **argv){

	cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║ \033[1m Time Series analysis for TETHYS                                       \033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";       	

	int N = atoi(argv[1]);	
	
	float S = atof(argv[3]);	
	
	float *Time;			
	Time =(float*) calloc (N,sizeof(float));
	
	float *in_den_0;			
	in_den_0 =(float*) calloc (N,sizeof(float));
	
	float *in_den_L;			
	in_den_L =(float*) calloc (N,sizeof(float));
	
	float *in_vel_0;			
	in_vel_0 =(float*) calloc (N,sizeof(float));
	
	float *in_vel_L;			
	in_vel_L =(float*) calloc (N,sizeof(float));

	//////////////////Reading input file ///////////////////////////////
	ifstream input;
	input.open(argv[2]);

	cout << "Reading input file";
	
	if(input.is_open())
	{
		int i=0;
		while(input.good())
		{	
			input >> Time[i] >> in_den_L[i] >> in_vel_L[i] >> in_den_0[i] >> in_vel_0[i];
			i++;	
		}
	}
	cout << "\tDONE!"<<"\tS="<<S<<endl;
	float dt = Time[2]- Time[1];
	cout << "dt " << dt <<endl;
	cout << "time max " << Time[N-1] <<endl;	
	////////////////////////////////////////////////////////////////////



	float saturation = 0.0;
	float tau = 0.0;
	float error = 0.0;
	ExtremaFinding(in_den_0, N, S, dt,saturation,tau, error, "den0.dat");
	cout << "Density saturation at x=0: " << saturation << endl;
	cout << "Time for 99% of saturation: " << tau <<endl;		

	ExtremaFinding(in_den_L, N, S, dt,saturation,tau, error, "denL.dat");
	cout << "Density saturation at x=L: " << saturation << endl;
	cout << "Time for 99% of saturation: " << tau <<endl;		

	ExtremaFinding(in_vel_0, N, S, dt,saturation,tau, error, "vel0.dat");
	cout << "Velocity saturation at x=0: " << saturation << endl
	cout << "Time for 99% of saturation: " << tau <<endl;		;

	ExtremaFinding(in_vel_L, N, S, dt,saturation,tau, error, "velL.dat");	
	cout << "Velocity saturation at x=L: " << saturation << endl;
	cout << "Time for 99% of saturation: " << tau <<endl;		
	
	
	
	
	
	return 0;
	
}
