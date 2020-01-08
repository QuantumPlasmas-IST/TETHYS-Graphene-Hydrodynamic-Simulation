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

#include "GraphHydro.h"
#include "TethysLib.h"


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

#ifndef C_SPEED
#    define C_SPEED 1000.0
#endif



using namespace std;




float PhaseVel(float sound, float fermi){
	float vel_phs = sqrt(sound*sound+0.5*fermi*fermi + 0.0625 );
	return vel_phs ;
}

float RealFreq(float sound, float fermi, float col_freq, int mode){
	float real;
	float vel_phs = PhaseVel(sound,fermi);
	float vel_phs_sqr = vel_phs*vel_phs ;
	if (1 < vel_phs ){
		mode = 2*mode-1;
	 }
	else{
		mode = 2*mode;
		}
	real =  fabs(vel_phs_sqr - 0.5625 ) * MAT_PI * mode / (2.0 * vel_phs );
	return real;
}
	
	
float ImagFreq(float sound, float fermi, float col_freq){
	float imag;	
	float vel_phs = PhaseVel(sound,fermi);
	float vel_phs_sqr = vel_phs*vel_phs ;
	imag = (vel_phs_sqr - 0.5625 ) * log(fabs( (vel_phs+0.75)/(vel_phs-0.75) )) / (2.0 * vel_phs ) - col_freq*(1-0.125/vel_phs);
	return imag;
}	

float TimeStepCFL(float dx, float sound, float fermi){
float dt;	

	if(fermi<10 && (sound-fermi) <= 3)
		dt = 0.5 * dx / (2*sound+sqrt(3*fermi*fermi + 24*sound*sound));
	else if (fermi<10 && (sound-fermi<= 10 - fermi))
		dt = 1.5 * dx / (2*sound+sqrt(3*fermi*fermi + 24*sound*sound));
	else if (fermi<15 && (sound-fermi<= 5))
		dt = 2 * dx / (2*sound+sqrt(3*fermi*fermi + 24*sound*sound));
	else if (fermi<30 && (sound-fermi<= 3))
		dt = 3 * dx / (2*sound+sqrt(3*fermi*fermi + 24*sound*sound));
	else
		dt = 4 * dx / (2*sound+sqrt(3*fermi*fermi + 24*sound*sound));
	
return dt;
}

void BoundaryCond(int type, int N, float * den, float * vel ){
	/*---------------*\
	| Free        | 1 |
	| Periodic 	  | 2 |
	| DS Boundary | 3 | 
	| DS+Driving  | 4 | 
	\*---------------*/
	switch(type){
		case 1 : den[0] = den[1];
				 den[N-1] = den[N-2];
				 vel[0] = vel[1];
				 vel[N-1] = vel[N-2];		
			break;
		case 2 : den[0] = den[N-2];
				 den[N-1] = den[1];
				 vel[0] = vel[N-2];
				 vel[N-1] = vel[1];	
			break;
		case 3 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[N-1] = den[N-2];
				 vel[N-1] = 1.0/den[N-1];			
			break;	
		case 4 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[N-1] = den[N-2];
				 vel[N-1] = (1.0 + 0.75*vel[1]*den[1])/den[N-1];			
			break;		
		default : den[0] = 1.0;
				  vel[0] = vel[1];
				  den[N-1] = den[N-2];
				  vel[N-1] = 1.0/den[N-1];			
	}
}


void InitialCondSine(int N, float dx,  float * den, float * vel){
  float L = (N-1)*dx;
  for (int i = 0; i < N; i++ )
  {
		float x = (float)i * dx;
		den[i] = 1.0 + 0.05*sin ( MAT_PI * x / L );
		vel[i] = 0.0;	
  }
}

void InitialCondRand(int N, float dx,  float * den, float * vel){
  srand (time(NULL));   
  for (int i = 0; i < N; i++ )
  {
		float noise = (float) rand()/ (float) RAND_MAX ;
		den[i] = 1.0 + 0.005*(noise-0.5);
		vel[i] = 0.0;
  }
}


void ExtremaFinding(float *vec_in, int N, float sound, float dt,float & sat, float & tau , float & error, std::string extremafile){		
	ofstream data_extrema;
	data_extrema.open(extremafile);
	
	data_extrema << "#pos_max" << "\t" << "Max" <<"\t"<< "pos_min" <<"\t"<< "Min"<<endl;	 
	
	int W = floor( 1.2*2*MAT_PI/(RealFreq(sound, 1.0, 1.0, 1)*dt));	
	int k = 0;
	int M = ceil(0.5*dt*N*RealFreq(sound, 1.0, 1.0, 1)/MAT_PI);
	float *vec_max;			
	vec_max =(float*) calloc (M,sizeof(float));
	float *vec_pos;			
	vec_pos =(float*) calloc (M,sizeof(float));
	
	for(int shift=0; shift < N-W ; shift += W){
		float maximum =  *max_element( vec_in + shift ,vec_in + shift + W );
		int pos_max =   max_element( vec_in + shift ,vec_in + shift + W ) - vec_in; 
		
		float minimum =  *min_element( vec_in + shift ,vec_in + shift + W );
		int pos_min =   min_element( vec_in + shift ,vec_in + shift + W ) - vec_in; 
		
		data_extrema << pos_max*dt << "\t" << maximum <<"\t"<< pos_min*dt <<"\t"<< minimum  <<endl;	 	
		vec_max[k] = maximum;
		vec_pos[k] = pos_max*dt;
		k++;
	}
	sat =  *max_element( vec_max ,vec_max + M );	
	for(int i=1;i<M-1;i++){
		if( vec_max[i] > 0.99*sat ){
			tau  = 	vec_pos[i];
			error = 0.5*(vec_pos[i+1]-vec_pos[i]);
			break;
		}
	}
	data_extrema.close();		
}

void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile){
	ofstream data_shock;
	data_shock.open(shockfile);
	for ( int i = 0; i < N; i++ )
	{
		if(i>=1){
			float schockD=0.0;
			schockD = (in[i]-in[i-1])/dx;
			if(abs(schockD) >=20){
				data_shock  <<t<<"\t"<< i*dx <<endl;	
			}
		}
	}
	data_shock.close();	
}
