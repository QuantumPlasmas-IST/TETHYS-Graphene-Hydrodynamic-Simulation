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

#include "Tethys1DLib.h"
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

GrapheneFluid1D::GrapheneFluid1D(int sizeN,float VELSND, float FERMI,float VISCO,float COL): Fluid1D(sizeN, VELSND, VISCO){
	vel_fer =FERMI;							
	col_freq =COL; 
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}

/*....................................................................*/	
/*.......... 1 Dimensional Fluid Class ...............................*/	
/*....................................................................*/	
Fluid1D::Fluid1D(int sizeN,float VELSND, float VISCO): TETHYSBase{sizeN}{	
	Nx = sizeN;
	vel_snd =VELSND;
	kin_vis =VISCO;
	den     = new float[sizeN]();
	vel     = new float[sizeN]();
	grad_vel= new float[sizeN]();
	cur     = new float[sizeN]();
	den_cor = new float[sizeN]();
	vel_cor = new float[sizeN]();
	cur_cor = new float[sizeN]();
	den_mid = new float[sizeN-1]();
	vel_mid = new float[sizeN-1]();
	grad_vel_mid = new float[sizeN-1]();
	vel_snd_arr = new float[sizeN-1]();
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
}	
	
Fluid1D::~Fluid1D(){			
	delete [] den;
	delete [] vel ;
	delete [] cur ;
	delete [] den_mid ;
	delete [] vel_mid ;
	delete [] den_cor ;
	delete [] vel_cor ;
	delete [] cur_cor ;
	delete [] vel_snd_arr ;
}

float  Fluid1D::DensityFlux(float n,float v,float S){
	float f1;
	f1 = n*v;
	return f1;		
}
float  Fluid1D::VelocityFlux(float n,float v,float dv,float S){
	float f2;
	f2 = 0.5*v*v + n - kin_vis*dv; 
	return f2;
}
float  Fluid1D::DensitySource(float n,float v,float S){
	return 0;
}
float  Fluid1D::VelocitySource(float n,float v,float S){
	return 0;
}	

void Fluid1D::CFLCondition(){
		dx = leng / ( float ) ( Nx - 1 );
		dt = dx/10.0;
}

void Fluid1D::SetSimulationTime(){
	Tmax=5+0.02*vel_snd+20.0/vel_snd;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]=SoundVelocityAnisotropy(i,dx,vel_snd);
	}
}
		
void Fluid1D::InitialCondRand(){
  srand (time(NULL));   
  for (int i = 0; i < Nx; i++ )
  {
		float noise = (float) rand()/ (float) RAND_MAX ;
		den[i] = 1.0 + 0.005*(noise-0.5);
  }	
}

void Fluid1D::SetKinVis(float x){ kin_vis=x;}
void Fluid1D::SetVelSnd(float x){ vel_snd=x; }
float Fluid1D::GetVelSnd(){ return vel_snd; }
float Fluid1D::GetKinVis(){ return kin_vis; }
float Fluid1D::GetTmax(){return Tmax;}
void Fluid1D::SetTmax(float x){ Tmax=x;}
float Fluid1D::GetDx(){return dx;}
float Fluid1D::GetDt(){return dt;}
int Fluid1D::SizeX(){ return Nx; }

void Fluid1D::Smooth(int width){
	AverageFilter( den ,den_cor, Nx, width);	
	AverageFilter( vel ,vel_cor, Nx, width);
	AverageFilter( cur ,cur_cor, Nx , width);
}



void Fluid1D::CreateFluidFile(){
	std::string previewfile = "preview_1D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid1D::WriteFluidFile(float t){
data_preview <<t<<"\t"<< den_cor[Nx-1] <<"\t"<< vel_cor[Nx-1] <<"\t"<< den_cor[0] <<"\t" << vel_cor[0] <<"\n";
}

void Fluid1D::Richtmyer(){
		//
		//  Calculating the velocity gradient 
		//
		for ( int i = 1; i < Nx-1 ; i++ )
		{
			grad_vel[i] = (-0.5*vel[i-1]+0.5*vel[i+1])/dx;
			if( i != Nx-1){
			grad_vel_mid[i] =(-0.5*vel_mid[i-1]+0.5*vel_mid[i+1])/dx;
			}
		}
		grad_vel[0] = (-1.5*vel[0]+2.0*vel[1]-0.5*vel[2])/dx;;
		grad_vel_mid[0] = (-1.5*vel_mid[0]+2.0*vel_mid[1]-0.5*vel_mid[2])/dx;
		grad_vel[Nx-1] =  ( 0.5*vel[Nx-1-2]-2.0*vel[Nx-1-1]+1.5*vel[Nx-1])/dx;
		grad_vel_mid[(Nx-1)-1] = ( 0.5*vel[(Nx-1)-1-2]-2.0*vel[(Nx-1)-1-1]+1.5*vel[(Nx-1)-1])/dx;
    	//
		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i < Nx - 1; i++ )
		{
			den_mid[i] = 0.5*( den[i] + den[i+1] )
				- ( 0.5*dt/dx ) * ( DensityFlux(den[i+1],vel[i+1],vel_snd_arr[i]) - DensityFlux(den[i],vel[i],vel_snd_arr[i]) ) 
				+ ( 0.5*dt    ) * DensitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),vel_snd_arr[i]) ;
			vel_mid[i] = 0.5*( vel[i] + vel[i+1] )
				- ( 0.5*dt/dx ) * ( VelocityFlux(den[i+1],vel[i+1],grad_vel[i+1],vel_snd_arr[i]) - VelocityFlux(den[i],vel[i],grad_vel[i],vel_snd_arr[i]) ) 
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),vel_snd_arr[i]) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],vel_snd_arr[i]) - DensityFlux(den_mid[i-1],vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * DensitySource(den[i],vel[i],vel_snd_arr[i]);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],grad_vel_mid[i],vel_snd_arr[i]) - VelocityFlux(den_mid[i-1],vel_mid[i-1],grad_vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * VelocitySource(den[i],vel[i],vel_snd_arr[i]);
			cur[i] = vel[i]*den[i];
		}
} 

/*....................................................................*/	
/*............ Derived Graphene Class  ...............................*/	
/*....................................................................*/	

float GrapheneFluid1D::DensityFlux(float n,float v,float S){
	float f1;
	f1 = n*v;
	return f1;			
}
float GrapheneFluid1D::VelocityFlux(float n,float v,float dv,float S){
	float f2;
	f2 = 0.25*v*v + vel_fer*vel_fer*0.5*log(n) + 2*S*S*sqrt(n) - kin_vis*dv; 
	return f2;			
}
float GrapheneFluid1D::DensitySource(float n,float v,float S){
	float Q1=0.0;
	return Q1;				
}
float GrapheneFluid1D::VelocitySource(float n,float v,float S){
	float Q2=0.0;
	Q2=-1.0*col_freq*(v-1);
	return Q2;			
}

void GrapheneFluid1D::SetVelFer(float x){ vel_fer=x; }
float GrapheneFluid1D::GetVelFer(){ return vel_fer; }
void GrapheneFluid1D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid1D::GetColFreq(){ return col_freq; }

void GrapheneFluid1D::WriteAtributes(){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	int total_steps=Tmax/dt;
	//Create the data space for the attribute.
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = grp_dat->createAttribute( "S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = grp_dat->createAttribute( "Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_kin_vis = grp_dat->createAttribute( "Kinetic viscosity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = grp_dat->createAttribute( "Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = grp_dat->createAttribute( "Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = grp_dat->createAttribute( "Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = grp_dat->createAttribute( "Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_space_points = grp_dat->createAttribute( "Number of spatial points", hdf5_int, atr_dataspace);
	Attribute atr_num_time_steps = grp_dat->createAttribute( "Number of time steps", hdf5_int, atr_dataspace);
	// Write the attribute data.
	atr_vel_snd.write( hdf5_float, &vel_snd);
	atr_vel_fer.write( hdf5_float, &vel_fer);
	atr_col_freq.write(hdf5_float, &col_freq);
	atr_kin_vis.write(hdf5_float, &kin_vis); 
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_num_space_points.write( hdf5_int, &Nx);
	atr_total_time.write( hdf5_float, &Tmax);
	atr_num_time_steps.write(hdf5_int, &total_steps);
	// Close the attributes.
	atr_num_time_steps.close();
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
}

void GrapheneFluid1D::CFLCondition(){
	dx = leng / ( float ) ( Nx - 1 );	
	float lambda;
	if(vel_snd<0.36*vel_fer){
		lambda=1.2*vel_fer;
	}else{
		lambda=1.97*vel_snd + 0.5*vel_fer;
	}
	dt = dx/lambda;				
}	



void GrapheneFluid1D::BoundaryCond(int type){
	/*---------------*\
	| Free        | 1 |
	| Periodic 	  | 2 |
	| DS Boundary | 3 | 
	| DS+Driving  | 4 | 
	\*---------------*/
	switch(type){
		case 1 : den[0] = den[1];
				 den[Nx-1] = den[Nx-2];
				 vel[0] = vel[1];
				 vel[Nx-1] = vel[Nx-2];		
			break;
		case 2 : den[0] = den[Nx-2];
				 den[Nx-1] = den[1];
				 vel[0] = vel[Nx-2];
				 vel[Nx-1] = vel[1];	
			break;
		case 3 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[Nx-1] = den[Nx-2];
				 vel[Nx-1] = 1.0/den[Nx-1];			
			break;	
		case 4 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[Nx-1] = den[Nx-2];
				 vel[Nx-1] = (1.0 + 0.75*vel[1]*den[1])/den[Nx-1];			
			break;		
		default : den[0] = 1.0;
				  vel[0] = vel[1];
				  den[Nx-1] = den[Nx-2];
				  vel[Nx-1] = 1.0/den[Nx-1];			
	}
}

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

void ElectroAnalysis::CreateElectroFile(GrapheneFluid1D& graphene){
	//graphene.SetFileName();
	std::string infix = graphene.GetInfix();
	std::string electrofile = "electro_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific; 	
}

void ElectroAnalysis::WriteElectroFile(float t,GrapheneFluid1D& graphene){
		float Q_net = this->NetCharge(graphene);
		float I_avg = this->AverageCurrent(graphene); 
		float P_ohm = this->OhmPower(graphene);
		float Dipole_var=this->ElectricDipoleVariation(graphene);
		float Dipole=this->ElectricDipole(graphene);
		data_electro <<t<<"\t"<< Q_net<<"\t"<<I_avg<<"\t"<<Q_net*Q_net*0.5 <<"\t"<<P_ohm<<"\t"<<Dipole<<"\t"<< Dipole_var <<"\n";
}

float ElectroAnalysis::NetCharge(GrapheneFluid1D& graphene){
	return Integral1D(graphene.SizeX(), graphene.GetDx(), graphene.den_cor);
}

float ElectroAnalysis::AverageCurrent(GrapheneFluid1D& graphene){
	return Integral1D(graphene.SizeX(), graphene.GetDx(), graphene.cur_cor);
}

float ElectroAnalysis::ElectricDipoleVariation(GrapheneFluid1D& graphene){
	return Integral1D(graphene.SizeX(), graphene.GetDx(), graphene.cur_cor);
}

float ElectroAnalysis::ElectricDipole(GrapheneFluid1D& graphene){
	float p=0.0;
	float dx=graphene.GetDx();
	for(int j=1;j<graphene.SizeX()/2;j++){	
		p += dx*(2*j-2)*graphene.den_cor[2*j-2] + 4*dx*(2*j-1)*graphene.den_cor[2*j-1] + dx*(2*j)*graphene.den_cor[2*j];
	}
	p = p*graphene.GetDx()/3.0;
	return p;
}

float ElectroAnalysis::OhmPower(GrapheneFluid1D& graphene){
	float itg=0.0;
	for(int j=1;j<graphene.SizeX()/2;j++){
		itg += graphene.cur_cor[2*j-2]*graphene.vel_cor[2*j-2] + 4*graphene.cur_cor[2*j-1]*graphene.vel_cor[2*j-1] + graphene.cur_cor[2*j]*graphene.vel_cor[2*j];
	}
	itg = itg*graphene.GetDx()/3.0;
	return itg;	
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
	grp_vel->close();
	hdf5file->close();
}

TETHYSBase::TETHYSBase(int sizeN){
	Nx = sizeN;
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);	
	char buffer [50];
	sprintf (buffer, "Fluido1D_Nx=%d", Nx);
	file_infix = buffer;
}




void TETHYSBase::CreateHDF5File(){

	std::string hdf5name = "hdf5_1D_" + this->GetInfix() + ".h5" ;	
	H5std_string  FILE_NAME( hdf5name );

    hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));
	grp_vel = new Group( hdf5file->createGroup( "/Data/Velocity" ));
	grp_cur = new Group( hdf5file->createGroup( "/Data/Current" ));
	
	
	/* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	//hsize_t     dimsf[2];              // dataset dimensions
	hsize_t     dimsf[1];              // dataset dimensions
	dimsf[0] = Nx;
	//dimsf[1] = Ny;
	
	dataspace_den = new DataSpace( RANK, dimsf );
	dataspace_vel = new DataSpace( RANK, dimsf );
	dataspace_cur = new DataSpace( RANK, dimsf );
	
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
