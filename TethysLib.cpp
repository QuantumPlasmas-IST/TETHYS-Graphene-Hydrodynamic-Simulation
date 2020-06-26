#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   
#include <cassert>

#include "TethysLib.h"
#include <H5Cpp.h>

using namespace H5;
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



float SoundVelocityAnisotropy(float i, float dx,float S){
	return S;
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





/*....................................................................*/
/*........ General Functions .........................................*/
/*....................................................................*/
void TETHYSBase::BannerDisplay(void){
cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║\033[2m  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆▆▆▆▆▆  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▖   ▗▆▆▆ ▗▆▆▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m  █▘  ▐█▌  ▝█  ▐█▌    ▝█  █▘  ▐█▌  ▝█  ▐█▌   ▐█▌     █▌ ▐█   ▐█▌     ▝█  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌▆▆▆█        ▐█▌      ▐█▌▆▆▆▐█▌      ▐█▌     ▝██▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌    ▗▉      ▐█▌      ▐█▌   ▐█▌      ▐█▌    ▗       ██  \033[0m║\n";
	cout<<"║\033[2m     ▆███▆    ▆███▆▆▆██▉     ▆███▆    ▆███▆ ▆███▆    ▆███▆   ▐█▆▆▆▆▆██▘  \033[0m║\n";
	cout<<"║                                                                         ║\n";
	cout<<"║ \033[1mTwo-dimensional Emitter of THz, Hydrodynamic Simulation.  Version 1.3.2\033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";                                                                                                                                                                                          
}

void TETHYSBase::WellcomeScreen(float vel_snd, float vel_fer, float col_freq,float viscosity, float dt,float dx, float Tmax){
	cout << "\nFermi velocity\t\033[1mvF\t"<< vel_fer <<" v\342\202\200\033[0m\n";
	if ( PhaseVel(vel_snd, vel_fer) < vel_fer){
		cout << "Phase velocity\t\033[1mS'\t" << PhaseVel(vel_snd, vel_fer)<<" v\342\202\200\033[0m  \033[1;5;7;31m WARNING plasmon in damping region \033[0m" <<endl;
	}else{
		cout << "Phase velocity\t\033[1mS'\t" << PhaseVel(vel_snd, vel_fer)<<" v\342\202\200\033[0m\n";
	}
	cout << "Viscosity \t\033[1m\316\267\t"<< viscosity <<"\033[0m\n";
	if (viscosity !=0.0){ cout << "Reynolds n. \t\033[1mRe\t"<< 1.0/viscosity <<"\033[0m\n";}
	cout << "Collision \t\033[1m\316\275\t"<< col_freq <<" v\342\202\200/L\n\033[0m\n";
	cout << "Theoretical frequency \033[1m\317\211=\317\211'+i\317\211''\033[0m\n";
	cout << "\033[1m\317\211'\t"<< RealFreq(vel_snd,vel_fer,col_freq,1) << " v\342\202\200/L\t2\317\200/\317\211'\t"<< 2.0*MAT_PI/RealFreq(vel_snd,vel_fer,col_freq,1)  << " L/v\342\202\200\033[0m\n";
	cout << "\033[1m\317\211''\t"<< ImagFreq(vel_snd,vel_fer,col_freq) <<" v\342\202\200/L\t2\317\200/\317\211''\t"<< 2.0*MAT_PI/ImagFreq(vel_snd,vel_fer,col_freq) <<" L/v\342\202\200\033[0m\n";
	cout <<"Determined maximum simulated time\t\033[1m\nT\342\202\230\342\202\220\342\202\223\t" <<Tmax<<" L/v\342\202\200\033[0m\n";
	cout <<"Discretisation\n";
	cout <<"\033[1m\316\224t\t"<<dt<<" L/v\342\202\200\t\316\224x\t"<<dx<<" L\033[0m\n"<<endl;
}

std::string TETHYSBase::GetInfix(){return file_infix;}


 TETHYSBase::~TETHYSBase(){
	grp_dat->close(); 
	grp_den->close(); 
	grp_velX->close();
	grp_velY->close();
	hdf5file->close();
}

TETHYSBase::TETHYSBase(int sizeNX,int sizeNY,int dimension){
	Nx = sizeNX;
	Ny = sizeNY;
	RANK=dimension;
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);	
	char buffer [50];
	sprintf (buffer, "Fluido1D_Nx=%d", Nx);
	file_infix = buffer;
}




void TETHYSBase::CreateHDF5File(){
	if(RANK==1){
		std::string hdf5name = "hdf5_1D_" + this->GetInfix() + ".h5" ;	
	}
	if(RANK==2){
		std::string hdf5name = "hdf5_2D_" + this->GetInfix() + ".h5" ;	
	}
	std::string hdf5name = "hdf5_1D_" + this->GetInfix() + ".h5" ;	
	H5std_string  FILE_NAME( hdf5name );
    hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

	grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));		
	grp_velX = new Group( hdf5file->createGroup( "/Data/VelocityX" ));
	grp_velY = new Group( hdf5file->createGroup( "/Data/VelocityY" ));

	if(RANK==1){		
		hsize_t     dimsf[1];              // dataset dimensions
		dimsf[0] = Nx;		
		dataspace_den = new DataSpace( RANK, dimsf );
		dataspace_velX = new DataSpace( RANK, dimsf );
	}
	if(RANK==2){		
		hsize_t     dimsf[2];
		dimsf[0] = Nx;
		dimsf[1] = Ny;
		dataspace_den = new DataSpace( RANK, dimsf );
		dataspace_velX = new DataSpace( RANK, dimsf );
		dataspace_velY = new DataSpace( RANK, dimsf );
	}		
}


//======================================================================
//======================================================================
//======================================================================

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

float Integral1D(int N, float ds, float * f){
	float itg=0.0;
	
	for(int j=1;j<N/2;j++){
		itg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	itg = itg*ds/3.0;
	return itg;	
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
			//out[i] = out[i] * size;
		}							
	}
}

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


