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


#include "TethysLib.h"
//#include "GraphHydro.h"

using namespace std;


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

#ifndef C_SPEED
#    define C_SPEED 1000.0
#endif


void BannerDisplay(void){
cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║\033[2m  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆▆▆▆▆▆  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▖   ▗▆▆▆ ▗▆▆▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m  █▘  ▐█▌  ▝█  ▐█▌    ▝█  █▘  ▐█▌  ▝█  ▐█▌   ▐█▌     █▌ ▐█   ▐█▌     ▝█  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌▆▆▆█        ▐█▌      ▐█▌▆▆▆▐█▌      ▐█▌     ▝██▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌    ▗▉      ▐█▌      ▐█▌   ▐█▌      ▐█▌    ▗       ██  \033[0m║\n";
	cout<<"║\033[2m     ▆███▆    ▆███▆▆▆██▉     ▆███▆    ▆███▆ ▆███▆    ▆███▆   ▐█▆▆▆▆▆██▘  \033[0m║\n";
	cout<<"║                                                                         ║\n";
	cout<<"║ \033[1mTwo-dimensional Emitter of THz, Hydrodynamic Simulation.  Version 1.2.5\033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";                                                                                                                                                                                          
}

void WellcomeScreen(float vel_snd, float vel_fer, float col_freq, float dt,float dx, float Tmax){
	cout << "\nFermi velocity\t\033[1mvF\t"<< vel_fer <<" v\342\202\200\033[0m\n";
	if ( PhaseVel(vel_snd, vel_fer) < vel_fer){
		cout << "Phase velocity\t\033[1mS'\t" << PhaseVel(vel_snd, vel_fer)<<" v\342\202\200\033[0m  \033[1;5;7;31m WARNING plasmon in damping region \033[0m" <<endl;
	}else{
		cout << "Phase velocity\t\033[1mS'\t" << PhaseVel(vel_snd, vel_fer)<<" v\342\202\200\033[0m\n";
	}
	cout << "Collision \t\033[1m\316\275\t"<< col_freq <<" v\342\202\200/L\n\033[0m\n";
	cout << "Theoretical frequency \033[1m\317\211=\317\211'+i\317\211''\033[0m\n";
	cout << "\033[1m\317\211'\t"<< RealFreq(vel_snd,vel_fer,col_freq,1) << " v\342\202\200/L\t2\317\200/\317\211'\t"<< 2.0*MAT_PI/RealFreq(vel_snd,vel_fer,col_freq,1)  << " L/v\342\202\200\033[0m\n";
	cout << "\033[1m\317\211''\t"<< ImagFreq(vel_snd,vel_fer,col_freq) <<" v\342\202\200/L\t2\317\200/\317\211''\t"<< 2.0*MAT_PI/ImagFreq(vel_snd,vel_fer,col_freq) <<" L/v\342\202\200\033[0m\n";
	cout <<"Determined maximum simulated time\t\033[1m\nT\342\202\230\342\202\220\342\202\223\t" <<Tmax<<" L/v\342\202\200\033[0m\n";
	cout <<"Discretisation\n";
	cout <<"\033[1m\316\224t\t"<<dt<<" L/v\342\202\200\t\316\224x\t"<<dx<<" L\033[0m\n"<<endl;
}



void RecordLogFile(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float Tmax){
	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	char buffer [80];
	strftime (buffer,80,"%F %H:%M:%S\n",time_info);
	logfile << "\n#Simulation @ " << buffer ;
	logfile << "#parameters:\n";
	logfile << "#vel_snd \t vel_fer \t col_freq  \t w' \t w'' \n";
	logfile << vel_snd <<"\t"<<vel_fer<< "\t"<< col_freq<<"\t"<< RealFreq(vel_snd,vel_fer,col_freq,1) <<"\t"<< ImagFreq(vel_snd,vel_fer,col_freq) <<"\n";
	logfile << "#discretisation:\n";
    logfile << "#dt\tdx\tTmax\ttime steps\tspace points\n";
	logfile << dt<<"\t"<<dx<<"\t"<<Tmax<<"\t"<< (int) Tmax/dt <<"\t"<< (int) 1/dx <<endl;
}


void Autocorrelation(float * out_gamma ,float * in , int crop, int size){
	int M = size - crop;
	float in_crop[M];
	for(int k =0;k < M;k++){
		in_crop[k] = in[k+crop];		
	}
	float sum;
	for(int lag=0; lag < M ;lag++){		
		for(int t=0; t < M ;t++){
			sum += in_crop[ t ] * in_crop[ (t + lag)%M];
		}
		out_gamma[lag] = sum;
		sum=0.0;
	}
}


void AverageFilter(float * vec_in, float * vec_out, int size , int width ){
	for ( int i = 0; i < size; i++ ){		
		if(i>=width &&i<=size-1-width){
			for(int k = i-width; k <= i+width;k++){			
				vec_out[i] += vec_in[k]; 
			}	
			vec_out[i] = vec_out[i]/(2.0*width+1.0);
			}
		else{
			vec_out[i] =	vec_in[i] ;
		}	
	}	
}


float RootMeanSquare(int N, float dt, float * f){
	float rms=0.0;
	
	for(int j=1;j<N/2;j++){
		rms += f[2*j-2]*f[2*j-2] + 4*f[2*j-1]*f[2*j-1] + f[2*j]*f[2*j];
	}
	rms = rms*dt/3.0;
	rms = sqrt( rms/(N*dt)  );
	return rms;	
}

float SignalAverage(int N, float dt, float * f){
	float avg=0.0;
	
	for(int j=1;j<N/2;j++){
		avg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	avg = avg*dt/3.0;
	avg = avg/(N*dt);
	return avg;	
}

float GaussKernel(int position , float t){
	float g;
	
	g = exp(-0.5*position*position/t);
	g = g/(sqrt(2*MAT_PI*t));
	
	return g;	
}


float GaussKernelDerivative(int position , float t){
	float g;
	
	g = -position*exp(-0.5*position*position/t)/t;
	g = g/(sqrt(2*MAT_PI*t));
	
	return g;	
}


void ConvolveGauss(int type, float M, float t, float * in, float * out, int size){
	if(type==0){	
		for(int i=0;i<size;i++){
			if(i>=M && i<size-M){
			for(int k=-M;k<=M;k++){
					out[i]  += in[i-k]*GaussKernel(k,t);
				}
			}
		}					
	}
	if(type==1){
		for(int i=0;i<size;i++){
			if(i>=M && i<size-M){
			for(int k=-M;k<=M;k++){
					out[i] += in[i-k]*GaussKernelDerivative(k,t);
				}
			}
			out[i] = out[i] * size;
		}							
	}
}



void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out ){
	//second order method
	//f[tempo][posicao]
	for(int i=1;i<size_rows-1;i++){
		for(int j=0;j<size_cols;j++){
			df_out[i][j] = (-0.5*f_in[i-1][j]+0.5*f_in[i+1][j])/dt;
		}
	}

	for(int j=0;j<size_cols;j++){
		df_out[0][j]           = (-1.5*f_in[0][j]+2.0*f_in[1][j]-0.5*f_in[2][j])/dt;
		df_out[size_rows-1][j] = ( 0.5*f_in[size_rows-1-2][j]-2.0*f_in[size_rows-1-1][j]+1.5*f_in[size_rows-1][j])/dt;
	}
}


     
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out ){
	//second order method
	//f[tempo][posicao]
	for(int i=0;i<size_rows;i++){
		for(int j=1;j<size_cols-1;j++){
			df_out[i][j] = (-0.5*f_in[i][j-1]+0.5*f_in[i][j+1])/dt;
		}
	}
	for(int i=0;i<size_rows;i++){
		df_out[i][0]           = (-1.5*f_in[i][0]+2.0*f_in[i][1]-0.5*f_in[i][2])/dt;
		df_out[i][size_cols-1] = ( 0.5*f_in[i][size_cols-1-2]-2.0*f_in[i][size_cols-1-1]+1.5*f_in[i][size_cols-1])/dt;
	}
}


float DensityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f1;
	f1 = den*vel;
	return f1;
}

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f2;
	f2 = 0.25*vel*vel + vel_fer*vel_fer*0.5*log(den) + 2*vel_snd*vel_snd*sqrt(den); 
	return f2;
}

/*
float EnergyFlux(float den,float vel,float vel_snd,float vel_fer){
	float f3;
	f3 = vel*pow(den,1.5);
 	return f3;
}
*/

float DensitySource(float den,float vel,float vel_snd,float vel_fer){
	float Q1=0.0;
return Q1;	
}

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq){
	float Q2=0.0;
	Q2=-1.0*col_freq*(vel-1);
return Q2;
}

/*
float EnergySource(float den,float den_der,float vel,float vel_snd,float vel_fer){
	float Q3=0.0;
	Q3=pow(vel_snd/vel_fer,2)*den*vel*den_der;
return Q3;
}
*/




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

