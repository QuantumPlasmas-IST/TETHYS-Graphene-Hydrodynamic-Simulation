// 2D version

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
#include "Tethys2DLib.h"
#include <H5Cpp.h>

using namespace H5;
using namespace std;


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

Fluid2D::Fluid2D(int sizeNx, int sizeNy, float VELSND, float VISCO) : TETHYSBase{sizeNx,sizeNy,2}{		
	Nx = sizeNx;
	Ny = sizeNy;
	vel_snd =VELSND;
	kin_vis =VISCO;
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	// main grid variables Nx*Ny
	den 		= new float[Nx*Ny](); 
	velX 		= new float[Nx*Ny](); 
	velY 		= new float[Nx*Ny](); 
	flxX 		= new float[Nx*Ny](); 
	flxY 		= new float[Nx*Ny](); 
	curX 		= new float[Nx*Ny](); 
	curY 		= new float[Nx*Ny](); 
	vel_snd_arr	= new float[Nx*Ny](); 
	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();  
	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();
	// 2nd Aux. Grid X (Nx-1)*(Ny)
	den_mid_x	= new float[(Nx-1)*Ny]();
	flxX_mid_x	= new float[(Nx-1)*Ny]();
	flxY_mid_x	= new float[(Nx-1)*Ny]();
	// 2nd Aux. Grid Y (Nx)*(Ny-1)
	den_mid_y	= new float[Nx*(Ny-1)](); 
	flxX_mid_y	= new float[Nx*(Ny-1)](); 
	flxY_mid_y	= new float[Nx*(Ny-1)](); 
}
	
Fluid2D::~Fluid2D(){
	delete [] den;
	delete [] velX;
	delete [] velY;
	delete [] flxX;
	delete [] flxY;
	delete [] curX;
	delete [] curY;
	delete [] den_mid_x;
	delete [] flxX_mid_x;
	delete [] flxY_mid_x;		
	delete [] den_mid_y;
	delete [] flxX_mid_y;
	delete [] flxY_mid_y;		
	delete [] den_mid;
	delete [] flxX_mid;
	delete [] flxY_mid;		
	delete [] vel_snd_arr;
}



void Fluid2D::SetSound(){ 
	for(int i = 0; i<Nx  ;i++){
		for(int j=0; j<Ny ; j++){
			//vel_snd_arr[i+j*Nx]=SoundVelocityAnisotropy(i,dx,j,dy, vel_snd);
			vel_snd_arr[i+j*Nx]=vel_snd;
		}
	}
}


float Fluid2D::GetVelSnd(){ return vel_snd; }
void Fluid2D::SetVelSnd(float x){ vel_snd=x; }
float Fluid2D::GetKinVis(){ return kin_vis; }
void Fluid2D::SetKinVis(float x){ kin_vis=x;}
float Fluid2D::GetDx(){return dx;}
void Fluid2D::SetDx(float x){ dx=x;} 
float Fluid2D::GetDy(){return dy;}
void Fluid2D::SetDy(float x){ dy=x;}
float Fluid2D::GetDt(){return dt;}
void Fluid2D::SetDt(float x){ dt=x;}

void Fluid2D::InitialCondRand(){
  	srand (time(NULL));   
  	for (int i = 0; i < Nx; i++ ){
  		for (int j=0; j<Ny; j++){
  			float noise = (float) rand()/ (float) RAND_MAX ;
			den[i+j*Nx] = 1.0 + 0.005*(noise-0.5);
  		}
  	}	
}

void Fluid2D::InitialCondTEST(){
  	for (int i = 0; i < Nx; i++ ){
  		for (int j=0; j<Ny; j++){
			float densi;
			if(i>=80&&i<=120&&j>=80&&j<=120){
			densi=0.2;	
			}
			else{
			densi=0.0;	
			}
			den[i+j*Nx] = 1.0 + densi;
			velX[i+j*Nx] = 0.1;
  		}
  	}	
}


void Fluid2D::MassFluxToVelocity(){
	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){
			velX[i+j*Nx]=flxX[i+j*Nx]/den[i+j*Nx];
			velY[i+j*Nx]=flxY[i+j*Nx]/den[i+j*Nx];
			
			
		    curX[i+j*Nx] =velX[i+j*Nx]*den[i+j*Nx];
			curY[i+j*Nx] =velY[i+j*Nx]*den[i+j*Nx];			
		}
	}
}

void Fluid2D::Richtmyer(){
		//
		//  Obtaining the 2nd auxiliary grids from averaging 2 neighbour cells in main grid
		//
		// mid_x
		for(int i=0; i<Nx-1; i++){
			for(int j=0; j<Ny; j++){
				den_mid_x[i+j*Nx] = 0.5*(den[i+j*Nx] + den[i+1+j*Nx]);
				flxX_mid_x[i+j*Nx] = 0.5*(flxX[i+j*Nx] + flxX[i+1+j*Nx]);
				flxY_mid_x[i+j*Nx] = 0.5*(flxY[i+j*Nx] + flxY[i+1+j*Nx]);
			}
		}
		// mid_y
		for(int i=0; i<Nx; i++){
			for(int j=0; j<Ny-1; j++){
				den_mid_y[i+j*Nx] = 0.5*(den[i+j*Nx] + den[i+(j+1)*Nx]);
				flxX_mid_y[i+j*Nx] = 0.5*(flxX[i+j*Nx] + flxX[i+(j+1)*Nx]);
				flxY_mid_y[i+j*Nx] = 0.5*(flxY[i+j*Nx] + flxY[i+(j+1)*Nx]);
			}
		}
	  	//
		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for(int i=0; i<Nx-1; i++){
			for(int j=0; j<Ny-1; j++){
				den_mid[i+j*Nx] = 0.25*(den[i+j*Nx] + den[i+1+j*Nx] + den[i+(j+1)*Nx] + den[i+1+(j+1)*Nx]) // How shall we include vel_snd_arr ?
								-0.5*dt*(
									DensityFluxX(den_mid_y[i+1+j*Nx], flxX_mid_y[i+1+j*Nx], flxY_mid_y[i+1+j*Nx],vel_snd_arr[i+j*Nx])-
									DensityFluxX(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx],vel_snd_arr[i+j*Nx]))/dx
								-0.5*dt*(
									DensityFluxY(den_mid_x[i+(j+1)*Nx], flxX_mid_x[i+(j+1)*Nx], flxY_mid_x[i+(j+1)*Nx],vel_snd_arr[i+j*Nx])-
									DensityFluxY(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx],vel_snd_arr[i+j*Nx]))/dy;
				flxX_mid[i+j*Nx] = 0.25*(flxX[i+j*Nx] + flxX[i+1+j*Nx] + flxX[i+(j+1)*Nx] + flxX[i+1+(j+1)*Nx])
								-0.5*dt*(
									MassFluxXFluxX(den_mid_y[i+1+j*Nx], flxX_mid_y[i+1+j*Nx], flxY_mid_y[i+1+j*Nx],vel_snd_arr[i+j*Nx])-
									MassFluxXFluxX(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx],vel_snd_arr[i+j*Nx]))/dx
								-0.5*dt*(
									MassFluxXFluxY(den_mid_x[i+(j+1)*Nx], flxX_mid_x[i+(j+1)*Nx], flxY_mid_x[i+(j+1)*Nx],vel_snd_arr[i+j*Nx])-
									MassFluxXFluxY(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx],vel_snd_arr[i+j*Nx]))/dy;
				flxY_mid[i+j*Nx] = 0.25*(flxY[i+j*Nx] + flxY[i+1+j*Nx] + flxY[i+(j+1)*Nx] + flxY[i+1+(j+1)*Nx])
								-0.5*dt*(
									MassFluxYFluxX(den_mid_y[i+1+j*Nx], flxX_mid_y[i+1+j*Nx], flxY_mid_y[i+1+j*Nx],vel_snd_arr[i+j*Nx])-
									MassFluxYFluxX(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx],vel_snd_arr[i+j*Nx]))/dx
								-0.5*dt*(
									MassFluxYFluxY(den_mid_x[i+(j+1)*Nx], flxX_mid_x[i+(j+1)*Nx], flxY_mid_x[i+(j+1)*Nx],vel_snd_arr[i+j*Nx])-
									MassFluxYFluxY(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx],vel_snd_arr[i+j*Nx]))/dy;
			}
		}
		//
		// Recalculate the 2nd auxiliary grids from averaging 2 neighbour cells in Aux grid / Note that dimension of the mid_x|mid_y is decreased by 2 in y|x
		//
		// mid_x
		for(int i=0; i<Nx-1; i++){
			for(int j=0; j<Ny-2; j++){
				den_mid_x[i+j*Nx] = 0.5*(den_mid[i+j*Nx] + den_mid[i+(j+1)*Nx]);
				flxX_mid_x[i+j*Nx] = 0.5*(flxX_mid[i+j*Nx] + flxX_mid[i+(j+1)*Nx]);
				flxY_mid_x[i+j*Nx] = 0.5*(flxY_mid[i+j*Nx] + flxY_mid[i+(j+1)*Nx]);
			}
		}
		// mid_y
		for(int i=0; i<Nx-2; i++){
			for(int j=0; j<Ny-1; j++){
				den_mid_y[i+j*Nx] = 0.5*(den_mid[i+j*Nx] + den_mid[i+1+j*Nx]);
				flxX_mid_y[i+j*Nx] = 0.5*(flxX_mid[i+j*Nx] + flxX_mid[i+1+j*Nx]);
				flxY_mid_y[i+j*Nx] = 0.5*(flxY_mid[i+j*Nx] + flxY_mid[i+1+j*Nx]);
			}
		}
		//
		// Obtain main grid at time k+1 - Remaining  cells to be definied by 
		//
		for(int i=0; i<Nx-2; i++){
			for(int j=0; j<Ny-2; j++){
				den[i+1+(j+1)*Nx] = den[i+1+(j+1)*Nx]
								-dt*(
									DensityFluxX(den_mid_x[i+1+j*Nx], flxX_mid_x[i+1+j*Nx], flxY_mid_x[i+1+j*Nx], vel_snd_arr[i+j*Nx])-
									DensityFluxX(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx], vel_snd_arr[i+j*Nx]))/dx
								-dt*(DensityFluxY(den_mid_y[i+(j+1)*Nx], flxX_mid_y[i+(j+1)*Nx], flxY_mid_y[i+(j+1)*Nx], vel_snd_arr[i+j*Nx])-
									DensityFluxY(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx], vel_snd_arr[i+j*Nx]))/dy;
				flxX[i+1+(j+1)*Nx] = flxX[i+1+(j+1)*Nx]
								-dt*(
									MassFluxXFluxX(den_mid_x[i+1+j*Nx], flxX_mid_x[i+1+j*Nx], flxY_mid_x[i+1+j*Nx], vel_snd_arr[i+j*Nx])-
									MassFluxXFluxX(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx], vel_snd_arr[i+j*Nx]))/dx
								-dt*(MassFluxXFluxY(den_mid_y[i+(j+1)*Nx], flxX_mid_y[i+(j+1)*Nx], flxY_mid_y[i+(j+1)*Nx], vel_snd_arr[i+j*Nx])-
									MassFluxXFluxY(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx], vel_snd_arr[i+j*Nx]))/dy;
				flxY[i+1+(j+1)*Nx] = flxY[i+1+(j+1)*Nx]
								-dt*(
									MassFluxYFluxX(den_mid_x[i+1+j*Nx], flxX_mid_x[i+1+j*Nx], flxY_mid_x[i+1+j*Nx], vel_snd_arr[i+j*Nx])-
									MassFluxYFluxX(den_mid_x[i+j*Nx], flxX_mid_x[i+j*Nx], flxY_mid_x[i+j*Nx], vel_snd_arr[i+j*Nx]))/dx
								-dt*(MassFluxYFluxY(den_mid_y[i+(j+1)*Nx], flxX_mid_y[i+(j+1)*Nx], flxY_mid_y[i+(j+1)*Nx], vel_snd_arr[i+j*Nx])-
									MassFluxYFluxY(den_mid_y[i+j*Nx], flxX_mid_y[i+j*Nx], flxY_mid_y[i+j*Nx], vel_snd_arr[i+j*Nx]))/dy;							
			}
		}
} 



void Fluid2D::CFLCondition(){ 
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0;
}



float  Fluid2D::DensityFluxX(float n,float flxX, float flxY, float S){ 
	float f1;
	f1 = flxX;
	return f1;		
}
float  Fluid2D::DensityFluxY(float n,float flxX, float flxY, float S){ 
	float f1;
	f1 = flxY;
	return f1;		
}
float  Fluid2D::DensitySource(float n,float velX, float velY, float S){
	float Q1 =0;
	return Q1;
}
float  Fluid2D::MassFluxXFluxX(float n,float flxX, float flxY, float S){
	float f2;
	f2 = flxX*flxX/n +n; 
	return f2;
}
float  Fluid2D::MassFluxXFluxY(float n,float flxX, float flxY, float S){
	float f2; 
	f2 = flxX*flxY/n;
	return f2;
}
float  Fluid2D::MassFluxYFluxX(float n,float flxX, float flxY, float S){
	float f3;
	f3 = flxX*flxY/n;
	return f3;
}
float  Fluid2D::MassFluxYFluxY(float n,float flxX, float flxY, float S){
	float f3;
	f3 = flxY*flxY/n + n;
	return f3;
}
float  Fluid2D::MassFluxXSource(float n,float flxX, float flxY, float S){
	float Q2 =0;
	return Q2;
}
float  Fluid2D::MassFluxYSource(float n,float flxX, float flxY, float S){
	float Q3 =0;
	return Q3;
}	

void Fluid2D::SetFileName(){
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
}

void Fluid2D::CreateFluidFile(){
	this->SetFileName();
	std::string previewfile = "preview_2D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid2D::WriteFluidFile(float t){
int j=Ny/2;
data_preview <<t<<"\t"<< den[Nx-1+j*Nx] <<"\t"<< velX[Nx-1+j*Nx] <<"\t"<< den[0+j*Nx] <<"\t" << velX[0+j*Nx] <<"\n";
}



void Fluid2D::SetSimulationTime(){
	Tmax=5+0.02*vel_snd+20.0/vel_snd;
}

/***************************************************** GrapheneFluid *****************************************************\
\*************************************************************************************************************************/


GrapheneFluid2D::GrapheneFluid2D(int sizeNx,int sizeNy,float VELSND, float FERMI,float VISCO,float COL): Fluid2D(sizeNx,sizeNy, VELSND, VISCO){
	vel_fer =FERMI;							
	col_freq =COL; 
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}


void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0+0.02*s+20.0/s);
}

void GrapheneFluid2D::MassFluxToVelocity(){
	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){
			velX[i+j*Nx]=flxX[i+j*Nx]*pow(den[i+j*Nx],-1.5);
			velY[i+j*Nx]=flxY[i+j*Nx]*pow(den[i+j*Nx],-1.5);
		    curX[i+j*Nx] = velX[i+j*Nx]*den[i+j*Nx];
			curY[i+j*Nx] = velY[i+j*Nx]*den[i+j*Nx];			
		}
	}
}


void GrapheneFluid2D::SetVelFer(float x){ vel_fer=x;	}
float GrapheneFluid2D::GetVelFer(){ return vel_fer;  }
void GrapheneFluid2D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid2D::GetColFreq(){ return col_freq; }


void GrapheneFluid2D::CFLCondition(){ // Eventual redefinition 
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );					

	//dt = 2.4/(vel_snd*sqrt(25.0/(dx*dx)+16.0/(dy*dy)));

	float lambda;
	
	if(vel_snd<0.36*vel_fer){
		lambda=1.2*vel_fer;
	}else{
		lambda=1.97*vel_snd + 0.5*vel_fer;
	}
		
	dt = dx/lambda;				


}	

//~ void GrapheneFluid2D::BoundaryCond(int type){
	//~ //       X        Y 
	//~ //1:       DS + open
	//~ //2:       DS + periodic
	//~ //3:     open + open
	//~ //4:     open + periodic
	//~ //5: periodic + open
	//~ //6: periodic + periodic
	//~ switch(type){
		//~ case 1 : 		
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=1.0;               			//constant density at x=0
				//~ den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				//~ flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			//~ }
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+1*Nx];
				//~ flxX[i+0*Nx] = flxX[i+1*Nx];
				//~ flxY[i+0*Nx] = flxY[i+1*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			//~ }	 	
			//~ break;
		//~ case 2 :
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=1.0;               			//constant density at x=0
				//~ den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				//~ flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			//~ }
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+1*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			//~ }
			//~ break;
		//~ case 3 :
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=den[1+j*Nx];                	
				//~ den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				//~ flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			//~ }		
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+1*Nx];
				//~ flxX[i+0*Nx] = flxX[i+1*Nx];
				//~ flxY[i+0*Nx] = flxY[i+1*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			//~ }	 			
			//~ break;
		//~ case 4 :
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=den[1+j*Nx];                	
				//~ den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				//~ flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			//~ }		
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+1*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			//~ }		
			//~ break;
		//~ case 5 :
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				//~ den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				//~ flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			//~ }		
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+1*Nx];
				//~ flxX[i+0*Nx] = flxX[i+1*Nx];
				//~ flxY[i+0*Nx] = flxY[i+1*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			//~ }	 	
			//~ break;
		//~ case 6 :
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				//~ den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				//~ flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			//~ }		
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+1*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			//~ }
			//~ break;
		//~ default : 			
			//~ for(int j=0;j<Ny;j++){	
				//~ den[0+j*Nx]=1.0;               			//constant density at x=0
				//~ den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				//~ flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				//~ flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				//~ flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				//~ flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			//~ }
			//~ for (int i=0; i<Nx; i++){
				//~ den[i+0*Nx] = den[i+1*Nx];
				//~ flxX[i+0*Nx] = flxX[i+1*Nx];
				//~ flxY[i+0*Nx] = flxY[i+1*Nx];
				//~ den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				//~ flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				//~ flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			//~ }	 	
	//~ }
//~ }


float  GrapheneFluid2D::DensityFluxX(float n,float flxX, float flxY, float S){ // Double Please Review this João 
	float f1;
	f1 = flxX/sqrt(n);
	return f1;		
}
float  GrapheneFluid2D::DensityFluxY(float n,float flxX, float flxY, float S){ 
	float f1;
	f1 = flxY/sqrt(n);
	return f1;		
}
float  GrapheneFluid2D::MassFluxXFluxX(float n,float flxX, float flxY, float S){
	float f2;
	f2 = flxX*flxX*pow(n,-1.5) +vel_fer*vel_fer*pow(n,1.5)/3.0+0.5*S*S*n*n; 
	return f2;
}
float  GrapheneFluid2D::MassFluxXFluxY(float n,float flxX, float flxY, float S){
	float f2; 
	f2 = flxX*flxY*pow(n,-1.5);
	return f2;
}
float  GrapheneFluid2D::MassFluxYFluxX(float n,float flxX, float flxY, float S){
	float f3;
	f3 = flxX*flxY*pow(n,-1.5);
	return f3;
}
float  GrapheneFluid2D::MassFluxYFluxY(float n,float flxX, float flxY, float S){
	float f3;
	f3 = flxY*flxY*pow(n,-1.5) + vel_fer*vel_fer*pow(n,1.5)/3.0+0.5*S*S*n*n;
	return f3;
}


// Pedro: para ja nao vamos incluir sources 
float  GrapheneFluid2D::DensitySource(float n,float flxX, float flxY, float S){
	float Q1 =0;
	return Q1;
}
float  GrapheneFluid2D::MassFluxXSource(float n,float flxX, float flxY, float S){
	float Q2 =0;
	return Q2;
}
float  GrapheneFluid2D::MassFluxYSource(float n,float flxX, float flxY, float S){
	float Q3 =0;
	return Q3;
}	

void GrapheneFluid2D::WriteAtributes(){
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
