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

#include "Tethys2DLib.h"

using namespace std;


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


Fluid2D::Fluid2D(int sizeNx, int sizeNy, float width) : lengY(width){		
	Nx = sizeNx;
	Ny = sizeNy;
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
			vel_snd_arr[i+j*Nx]=SoundVelocityAnisotropy(i,dx,j,dy, vel_snd);
		}
	}
}


float Fluid2D::GetVelSnd(){ return vel_snd; }
void Fluid2D::SetVelSnd(float x){ vel_snd=x; }
float Fluid2D::GetKinVis(){ return kin_vis; }
void Fluid2D::SetKinVis(float x){ kin_vis=x;}
void Fluid2D::SetTmax(float x){ Tmax=x;}
float Fluid2D::GetTmax(){return Tmax;}
float Fluid2D::GetDx(){return dx;}
void Fluid2D::SetDx(float x){ dx=x;} 
float Fluid2D::GetDy(){return dy;}
void Fluid2D::SetDy(float x){ dy=x;}
float Fluid2D::GetDt(){return dt;}
void Fluid2D::SetDt(float x){ dt=x;}
int Fluid2D::SizeX(){ return Nx; }
int Fluid2D::SizeY(){ return Ny; }

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


/*void Fluid2D::Smooth(int width){
	AverageFilter( den ,den_cor, Nx, Ny, width);	
	AverageFilter( velX ,velX_cor, Nx, Ny, width);
	AverageFilter( curX ,curX_cor, Nx , Ny, width);
	AverageFilter( velY ,velY_cor, Ny, Ny, width);
	AverageFilter( curY ,curY_cor, Ny , Ny, width);
}*/

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

void Fluid2D::BoundaryCond(int type){	
	// For now, periodic end in y 
	for (int i=0; i<Nx; i++)
	{
		den[i+0*Nx] = den[i+(Ny-2)*Nx];
		flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
		flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
		
		den[i+(Ny-1)*Nx] = den[i+1*Nx];
		flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
		flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
	}
	// in x
	/*---------------*\
	| Free        | 1 |
	| Periodic 	  | 2 |
	\*---------------*/
	switch(type){
		case 1 : 
		for(int j=0; j<Ny-1; j++){
			den[0+j*Nx] = den[1+j*Nx];
			flxX[0+j*Nx] = flxX[1+j*Nx];
			flxY[0+j*Nx] = flxY[1+j*Nx];
			den[Nx-1+j*Nx] = den[Nx-2+j*Nx];
			flxX[Nx-1+j*Nx] = flxX[Nx-2+j*Nx];
			flxY[Nx-1+j*Nx] = flxY[Nx-2+j*Nx];	
		}	
			break;
		case 2 :
		for(int j=0; j<Ny-1; j++){
			den[0+j*Nx] = den[Nx-2+j*Nx];
			flxX[0+j*Nx] = flxX[Nx-2+j*Nx];
			flxY[0+j*Nx] = flxY[Nx-2+j*Nx];
			den[Nx-1+j*Nx] = den[1+j*Nx];
			flxX[Nx-1+j*Nx] = flxX[1+j*Nx];
			flxY[Nx-1+j*Nx] = flxY[1+j*Nx];	
		}
			break;	
		default : 
		for(int j=0; j<Ny-1; j++){
			den[0+j*Nx] = den[Nx-2+j*Nx];
			flxX[0+j*Nx] = flxX[Nx-2+j*Nx];
			flxY[0+j*Nx] = flxY[Nx-2+j*Nx];
			den[Nx-1+j*Nx] = den[1+j*Nx];
			flxX[Nx-1+j*Nx] = flxX[1+j*Nx];
			flxY[Nx-1+j*Nx] = flxY[1+j*Nx];	
		}		
	}
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


/***************************************************** GrapheneFluid *****************************************************\
\*************************************************************************************************************************/

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
	if(vel_fer<10 && (vel_snd-vel_fer) <= 3)
		dt = 0.5 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<10 && (vel_snd-vel_fer<= 10 - vel_fer))
		dt = 1.5 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<15 && (vel_snd-vel_fer<= 5))
		dt = 2 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<30 && (vel_snd-vel_fer<= 3))
		dt = 3 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else
		dt = 4 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
}	


void GrapheneFluid2D::BoundaryCond(int type){
	//       X        Y 
	//1:       DS + open
	//2:       DS + periodic
	//3:     open + open
	//4:     open + periodic
	//5: periodic + open
	//6: periodic + periodic
	switch(type){
		case 1 : 		
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
			break;
		case 2 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}
			break;
		case 3 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[1+j*Nx];                	
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 			
			break;
		case 4 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[1+j*Nx];                	
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}		
			break;
		case 5 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
			break;
		case 6 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}
			break;
		default : 			
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
	}
}


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


/**/

/**************************************************** Other Functions ****************************************************\
\*************************************************************************************************************************/



float SoundVelocityAnisotropy(float i, float dx, float j, float dy, float S){
	return S;
}


void AverageFilter(float * vec_in, float * vec_out, int sizeX, int sizeY, int width ){ 
			for ( int i = 0; i < sizeX; i++ ){
				for(int j = 0; j < sizeY; j++){
					vec_out[i+j*sizeX] = 0;
					if(i>=width && i<=sizeX-1-width && j>=width && j<=sizeY-1-width){
						for (int m = i-width; m <= i+width ; m++){
							for (int n = j-width; m <= j+width ; m++){
								vec_out[i+j*sizeX] += vec_in[m+n*sizeX];
							}
						}
						vec_out[i+j*sizeX]/=(2.0*width+1.0);
						vec_out[i+j*sizeX]/=(2.0*width+1.0);
					}
					else{
						vec_out[i+j*sizeX] =	vec_in[i+j*sizeX];
					}
				}	
			}	
		}




void BannerDisplay(void){
cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║\033[2m  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆▆▆▆▆▆  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▖   ▗▆▆▆ ▗▆▆▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m  █▘  ▐█▌  ▝█  ▐█▌    ▝█  █▘  ▐█▌  ▝█  ▐█▌   ▐█▌     █▌ ▐█   ▐█▌     ▝█  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌▆▆▆█        ▐█▌      ▐█▌▆▆▆▐█▌      ▐█▌     ▝██▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌    ▗▉      ▐█▌      ▐█▌   ▐█▌      ▐█▌    ▗       ██  \033[0m║\n";
	cout<<"║\033[2m     ▆███▆    ▆███▆▆▆██▉     ▆███▆    ▆███▆ ▆███▆    ▆███▆   ▐█▆▆▆▆▆██▘  \033[0m║\n";
	cout<<"║                                                                         ║\n";
	cout<<"║ \033[1mTwo-dimensional Emitter of THz, Hydrodynamic Simulation.  Version 2.0.0\033[0m ║\n";
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

