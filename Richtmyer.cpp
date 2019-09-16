#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   

#include "dyakonovshur.h"


using namespace std;


float F1(float n,float v,float s);

float F2(float n,float v,float s);
float F2sqrt(float n,float v,float s,float vf);


int main(int argc, char **argv){

	cout << "*******************************************************"<<endl;
	cout << "************** ELECTRON FLOW SIMULATION ***************"<<endl;
	cout << "**************     Richtmyer method     ***************"<<endl;
	cout << "*******************************************************"<<endl;


	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	logfile << "\n#Simulation @ " << asctime(timeinfo) ;
	
	
	int N=201; 								// number of spatial points
	float t=0.0,L=1.0;						// spatial Length
	float dx;								// spatial discretization
	float dt;								// time step
	float S;							    // Sound speed

	float *n;							
	n =(float*) calloc (N,sizeof(float));
	float *n_mid;							
	n_mid = (float*) calloc (N-1,sizeof(float));
	float *v;							
 	v = (float*) calloc (N,sizeof(float));
	float *v_mid;						
	v_mid = (float*) calloc (N-1,sizeof(float));
	
	float *n_cor;							
	n_cor = (float*) calloc (N,sizeof(float));
	float *v_cor;						
	v_cor = (float*) calloc (N,sizeof(float));
 	float *j_cor;						
	j_cor = (float*) calloc (N,sizeof(float));
 	
 	
 	int flag=0;
	
 	
	if(argc!=1){
		
		S = atof(argv[1]);
			
		flag = atoi(argv[2]);	
	
		}
	else{
		cout << "Define S value: ";
		cin >> S;
		}
	
	
		
	dx = L / ( float ) ( N - 1 );
	
	
	
	if(S<5){
		dt = dx / (5*S);
	}
	else{
		if(S>8 && S<10){
			dt = dx / (30+3*S);
		}
		else{
			dt = dx / (20+2*S);		
			//dt = dx / (5+1.5*S);		
		}
	}
	
	float *s;							
	s =(float*) calloc (N,sizeof(float));	
	for(int i = 0; i<N  ;i++){
		//s[i]= S - 0.15*S*( dx*i- floor(dx*i) );
		s[i]=S;
	}

		string densityfile = "density_" + to_string(S)+ ".dat" ;
		densityfile.erase (densityfile.end()-9, densityfile.end()-5);
		
		string velocityfile = "velocity_" + to_string(S)+ ".dat" ;
		velocityfile.erase (velocityfile.end()-9, velocityfile.end()-5);
		
		string currentfile = "current_" + to_string(S)+ ".dat" ;
		currentfile.erase (currentfile.end()-9, currentfile.end()-5);
							
		ofstream data_density;
		data_density.open (densityfile);
		data_density << fixed ;
		data_density << setprecision(6);

		ofstream data_velocity;
		data_velocity.open (velocityfile);
		data_velocity << fixed ;
		data_velocity << setprecision(6);

		ofstream data_current;
		data_current.open (currentfile);
		data_current << fixed ;
		data_current << setprecision(6);	
	

	string electrofile = "electro_" + to_string(S)+ ".dat" ;
	electrofile.erase (electrofile.end()-9, electrofile.end()-5);

	string slicefile = "slice_" + to_string(S)+ ".dat" ;
	slicefile.erase (slicefile.end()-9, slicefile.end()-5);
	
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific; 

	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	

	
	
//	cout << "\n*******************************************************"<< endl;
	cout << "Sound speed S\t"<< S <<endl;
	cout <<"dt= "<<dt<<"\tdx= "<<dx<<endl;
	cout << "Predicted w'= "<< RealFreq(S,1.0,1.0,1) << "\t1/w'= "<< 1.0/RealFreq(S,1.0,1.0,1)  << endl;
	cout << "Predicted w''= "<< ImagFreq(S,1.0,1.0) <<"\t1/w''= "<< 1.0/ImagFreq(S,1.0,1.0) <<endl;
	
	logfile << "#S \t dt \t dx \t w' \t w'' " << endl;
	logfile << S <<"\t"<< dt <<"\t"<< dx <<"\t"<< RealFreq(S,1.0,1.0,1) <<"\t"<< ImagFreq(S,1.0,1.0) ;
	
//	cout << "*******************************************************"<< endl;
	
	float T_max=10.0;
	
	cout <<"Determined maximum simulated time\t" <<T_max<<endl;
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	//InitialCondRand(N, dx, n_cor, v_cor);
	//BoundaryCond(3, N, n_cor, v_cor);
	
	InitialCondRand(N, dx, n, v);
	BoundaryCond(3, N, n, v);
	////////////////////////////////////////////////////////////////////
	
	if(flag){
		for(int i = 0; i<N  ;i++)
		{
			data_density   <<  n[i] <<"\t";
			data_current   <<  v[i]*n[i] <<"\t";
			data_velocity  <<  v[i] <<"\t";
		}
	}
	
	cout << "Running"<<endl;
	
	int passo=0;
	
	while(t<=T_max && isfinite(v[(N-1)/2]))
	{	
		++passo;
		t += dt;
		
		if(flag && passo % 35 == 0 ){
			data_density  << "\n";
			data_current  << "\n";
			data_velocity << "\n";
		}
		
		//
		//  Half step calculate n and p at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i < N - 1; i++ )
		{
			n_mid[i] = 0.5*( n[i] + n[i+1] )
				- ( 0.5*dt/dx ) * ( F1(n[i+1],v[i+1],s[i+1]) - F1(n[i],v[i],s[i]) ) ;
			v_mid[i] = 0.5*( v[i] + v[i+1] )
				- ( 0.5*dt/dx ) * ( F2(n[i+1],v[i+1],s[i+1]) - F2(n[i],v[i],s[i]) ) ;	
		//	v_mid[i] = 0.5*( v[i] + v[i+1] )
		//		- ( 0.5*dt/dx ) * ( F2sqrt(n[i+1],v[i+1],s[i+1],5.0) - F2sqrt(n[i],v[i],s[i],5.0) ) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < N - 1; i++ )
		{
			n[i] = n[i] - (dt/dx) * ( F1(n_mid[i],v_mid[i],s[i]) - F1(n_mid[i-1],v_mid[i-1],s[i-1]) );
			v[i] = v[i] - (dt/dx) * ( F2(n_mid[i],v_mid[i],s[i]) - F2(n_mid[i-1],v_mid[i-1],s[i-1]) );
			//v[i] = v[i] - (dt/dx) * ( F2sqrt(n_mid[i],v_mid[i],s[i],5.0) - F2sqrt(n_mid[i-1],v_mid[i-1],s[i-1],5.0) );
		}
		
		BoundaryCond(3, N, n, v);
		
		AverageFilter( n ,n_cor, N , 2);	
		AverageFilter( v ,v_cor, N , 2);
		

			
			
			for ( int i = 0; i < N; i++ )
			{	
				
				
				j_cor[i] = 	v_cor[i]*n_cor[i];	
				if(flag){
				//Record full data
					if(passo % 35 == 0){
						data_density   << n_cor[i] <<"\t";
						data_current   << j_cor[i] <<"\t";
						data_velocity  << v_cor[i] <<"\t";
					}
				}
			}
		
		
		//Record end points
		data_slice <<t<<"\t"<< n_cor[N-1] <<"\t"<< v_cor[N-1] <<"\t"<< n_cor[0] <<"\t" << v_cor[0] <<endl;
		data_electro <<t<<"\t"<< n_cor[N-1]-1.0 <<"\t"<< n_cor[0]*v_cor[0] <<"\t"<<  TotalElectricDipole(N,dx,n_cor)<<"\t"<<  DtElectricDipole(N,dx,j_cor) <<"\t"<< KineticEnergy(N,dx, n_cor, v_cor)  <<endl;
	}
	cout << "DONE!" <<endl;
	cout << "*******************************************************"<<endl;

	free(n);
	free(n_mid);
	free(n_cor);
	free(v);
	free(v_mid);
	free(v_cor);
	free(j_cor);


		data_density.close();
		data_velocity.close();
		data_current.close();

	data_slice.close();
	data_electro.close();
	
	
	
	
	return 0;
}

float F1(float n,float v,float s){
	float f1;
	
	f1 = n*v;
	
	return f1;
}

float F2(float n,float v,float s){
	float f2;
	
	f2 = 0.5*v*v + s*s*n;
	
	return f2;
}

float F2sqrt(float n,float v,float s,float vf){
	float f2;
	
	f2 = 0.5*v*v + s*s*n +vf*vf*sqrt(n);
	
	return f2;
}

