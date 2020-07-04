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
Fluid1D::Fluid1D(int sizeNX,float VELSND, float VISCO): TETHYSBase{sizeNX,0,1}{	
	Nx = sizeNX;
	vel_snd =VELSND;
	kin_vis =VISCO;
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	den     = new float[sizeNX]();
	vel     = new float[sizeNX]();
	grad_vel= new float[sizeNX]();
	cur     = new float[sizeNX]();
	den_cor = new float[sizeNX]();
	vel_cor = new float[sizeNX]();
	cur_cor = new float[sizeNX]();
	den_mid = new float[sizeNX-1]();
	vel_mid = new float[sizeNX-1]();
	grad_vel_mid = new float[sizeNX-1]();
	vel_snd_arr = new float[sizeNX-1]();
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
	for(int i = 0; i<Nx-1  ;i++){
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
float Fluid1D::GetDx(){return dx;}
float Fluid1D::GetDt(){return dt;}


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
		//  Calculating the velocity gradient at k time
		//
		for ( int i = 1; i <= Nx-2 ; i++ )
		{
			grad_vel[i] = (-0.5*vel[i-1]+0.5*vel[i+1])/dx;
		}
		grad_vel[0] = (-1.5*vel[0]+2.0*vel[1]-0.5*vel[2])/dx;;
		grad_vel[Nx-1] =  ( 0.5*vel[Nx-1-2]-2.0*vel[Nx-1-1]+1.5*vel[Nx-1])/dx;

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
        //  Calculating the velocity gradient at k+1/2 time
        //
        for ( int i = 1; i <= Nx-3 ; i++ )
        {
            grad_vel_mid[i] =(-0.5*vel_mid[i-1]+0.5*vel_mid[i+1])/dx;
        }
        grad_vel_mid[0] = (-1.5*vel_mid[0]+2.0*vel_mid[1]-0.5*vel_mid[2])/dx;
        grad_vel_mid[(Nx-1)-1] = ( 0.5*vel[(Nx-1)-3]-2.0*vel[(Nx-1)-2]+1.5*vel[(Nx-1)-1])/dx;
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			float den_old = den[i];
			float vel_old = vel[i];
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],vel_snd_arr[i]) - DensityFlux(den_mid[i-1],vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * DensitySource(den_old,vel_old,vel_snd_arr[i]);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],grad_vel_mid[i],vel_snd_arr[i]) - VelocityFlux(den_mid[i-1],vel_mid[i-1],grad_vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * VelocitySource(den_old,vel_old,vel_snd_arr[i]);
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

