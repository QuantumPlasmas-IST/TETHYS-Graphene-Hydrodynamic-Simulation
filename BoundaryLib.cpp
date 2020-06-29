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
#include "Tethys2DLib.h"
#include "BoundaryLib.h"
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


//~ BoundaryCondition::BoundaryCondition(int sizeNx,int sizeNy,int dimensions)  : TETHYSBase{sizeNx,sizeNy,dimensions}{
//~ }




void BoundaryCondition::XFree(GrapheneFluid1D& graphene){
	int Nx=graphene.SizeX();
	graphene.den[0] = graphene.den[1];
	graphene.den[Nx-1] =  graphene.den[Nx-2];
	graphene.vel[0] = graphene.vel[1];
	graphene.vel[Nx-1] = graphene.vel[Nx-2];	
}
void BoundaryCondition::XPeriodic(GrapheneFluid1D& graphene){
	int Nx=graphene.SizeX();
	graphene.den[0] = graphene.den[Nx-2];
	graphene.den[Nx-1] = graphene.den[1];
	graphene.vel[0] = graphene.vel[Nx-2];
	graphene.vel[Nx-1] = graphene.vel[1];	
}
void BoundaryCondition::XFree(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	
	for(int j=0;j<Ny;j++){	
		graphene.den[0+j*Nx]=graphene.den[1+j*Nx];                	
		graphene.den[Nx-1+j*Nx]=graphene.den[Nx-2+j*Nx]; 			//free density at x=L
		graphene.flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
		graphene.flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
		graphene.flxX[0+j*Nx] = graphene.flxX[1+j*Nx]*pow(graphene.den[1+j*Nx],-1.5);			//free flux at x=0
		graphene.flxX[Nx-1+j*Nx] =  graphene.flxX[Nx-2+j*Nx];
	}	
}
void BoundaryCondition::XPeriodic(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	
	for(int j=0;j<Ny;j++){	
		graphene.den[0+j*Nx]=graphene.den[(Nx-2)+j*Nx];                	
		graphene.den[Nx-1+j*Nx]=graphene.den[1+j*Nx]; 			
		graphene.flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
		graphene.flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
		graphene.flxX[0+j*Nx] = graphene.flxX[(Ny-2)+j*Nx];
		graphene.flxX[Nx-1+j*Nx] =  graphene.flxX[1+j*Nx];
	}	
}
void BoundaryCondition::YFree(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	
	for (int i=0; i<Nx; i++){
		graphene.den[i+0*Nx] = graphene.den[i+1*Nx];
		graphene.flxX[i+0*Nx] = graphene.flxX[i+1*Nx];
		graphene.flxY[i+0*Nx] = graphene.flxY[i+1*Nx];
		graphene.den[i+(Ny-1)*Nx] = graphene.den[i+(Ny-2)*Nx];
		graphene.flxX[i+(Ny-1)*Nx] = graphene.flxX[i+(Ny-2)*Nx];
		graphene.flxY[i+(Ny-1)*Nx] = graphene.flxY[i+(Ny-2)*Nx];
	}	 	
}
void BoundaryCondition::YPeriodic(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();	
	
	for (int i=0; i<Nx; i++){
		graphene.den[i+0*Nx] = graphene.den[i+(Ny-2)*Nx];
		graphene.flxX[i+0*Nx] = graphene.flxX[i+(Ny-2)*Nx];
		graphene.flxY[i+0*Nx] = graphene.flxY[i+(Ny-2)*Nx];
		graphene.den[i+(Ny-1)*Nx] = graphene.den[i+1*Nx];
		graphene.flxX[i+(Ny-1)*Nx] = graphene.flxX[i+1*Nx];
		graphene.flxY[i+(Ny-1)*Nx] = graphene.flxY[i+1*Nx];
	}
}



void BoundaryCondition::Dirichlet::Density(GrapheneFluid1D& graphene, float L, float R){
	int Nx=graphene.SizeX();
	graphene.den[0] = L;
	graphene.den[Nx-1] = R;
}
void BoundaryCondition::Dirichlet::VelocityX(GrapheneFluid1D& graphene, float L, float R){
	int Nx=graphene.SizeX();
	graphene.vel[0] = L;
	graphene.vel[Nx-1] = R;
}  
void BoundaryCondition::Dirichlet::Density(GrapheneFluid2D& graphene, float L, float R, float T, float B){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int j=0; j<Ny; j++){
		graphene.den[0+j*Nx] = L;
		graphene.den[Nx-1+j*Nx] = R;
	}
	for (int i=0; i<Nx; i++){
		graphene.den[i+(Ny-1)*Nx] = T;		
		graphene.den[i+0*Nx] = B;
	}
} 

void BoundaryCondition::Dirichlet::VelocityX(GrapheneFluid2D& graphene, float L, float R, float T, float B){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int j=0; j<Ny; j++){
		graphene.velX[0+j*Nx] = L;
		graphene.velX[Nx-1+j*Nx] = R;
	}
	for (int i=0; i<Nx; i++){
		graphene.velX[i+(Ny-1)*Nx] = T;		
		graphene.velX[i+0*Nx] = B;
	}
} 
void BoundaryCondition::Dirichlet::VelocityY(GrapheneFluid2D& graphene, float L, float R, float T, float B){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int j=0; j<Ny; j++){
		graphene.velY[0+j*Nx] = L;
		graphene.velY[Nx-1+j*Nx] = R;
	}
	for (int i=0; i<Nx; i++){
		graphene.velY[i+(Ny-1)*Nx] = T;		
		graphene.velY[i+0*Nx] = B;
	}
} 

void BoundaryCondition::DyakonovShur::X(GrapheneFluid1D& graphene) {
	int Nx=graphene.SizeX();
	graphene.den[0] = 1.0;
	graphene.vel[0] = graphene.vel[1];
	graphene.den[Nx-1] = graphene.den[Nx-2];
	graphene.vel[Nx-1] = 1.0/graphene.den[Nx-1];		
}

void BoundaryCondition::DyakonovShur::X(GrapheneFluid2D& graphene) {
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	 
	for(int j=0;j<Ny;j++){	
		graphene.den[0+j*Nx]=1.0;               			//constant density at x=0
		graphene.den[Nx-1+j*Nx]=graphene.den[Nx-2+j*Nx]; 			//free density at x=L
		graphene.flxX[0+j*Nx] = graphene.flxX[1+j*Nx]*pow(graphene.den[1+j*Nx],-1.5);			//free flux at x=0
		graphene.flxX[Nx-1+j*Nx] = sqrt(graphene.den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
		graphene.flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
		graphene.flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
	}	
}

void BoundaryCondition::DyakonovShur::YFree(GrapheneFluid2D& graphene) {
	BoundaryCondition enclosing;	
	enclosing.YFree(graphene);
}

void BoundaryCondition::DyakonovShur::YPeriodic(GrapheneFluid2D& graphene) {
	BoundaryCondition enclosing;
	enclosing.YPeriodic(graphene);
}
