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
		int Left=0+j*Nx;
		int Right=Nx-1+j*Nx;
		graphene.den[Left]=graphene.den[Left+1];
		graphene.den[Right]=graphene.den[Right-1];			//free density at x=L
		graphene.flxY[Left] = 0.0f; 					//flux only on x at x=0
		graphene.flxY[Right] = 0.0f ;					//idem at x=L
		graphene.flxX[Left] = graphene.flxX[Left+1]*pow(graphene.den[Left+1],-1.5f);			//free flux at x=0
		graphene.flxX[Right] =  graphene.flxX[Right-1];
	}	
}
void BoundaryCondition::XPeriodic(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	
	for(int j=0;j<Ny;j++){
		int Left=0+j*Nx;
		int Right=Nx-1+j*Nx;
		graphene.den[Left]=graphene.den[Right-1];
		graphene.den[Right]=graphene.den[1+j*Nx];
		graphene.flxY[Left] = 0.0; 					//flux only on x at x=0
		graphene.flxY[Right] = 0.0 ;					//idem at x=L
		graphene.flxX[Left] = graphene.flxX[Right-1];
		graphene.flxX[Right] =  graphene.flxX[Left+1];
	}	
}
void BoundaryCondition::YFree(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int i=0; i<Nx; i++){
		int Bottom=i; //i+0*Nx
		int Top=i+(Ny-1)*Nx;
		graphene.den[Bottom] = graphene.den[Bottom+Nx];
		graphene.flxX[Bottom] = graphene.flxX[Bottom+Nx];
		graphene.flxY[Bottom] = graphene.flxY[Bottom+Nx];
		graphene.den[Top] = graphene.den[Top-Nx];
		graphene.flxX[Top] = graphene.flxX[Top-Nx];
		graphene.flxY[Top] = graphene.flxY[Top-Nx];
	}	 	
}
void BoundaryCondition::YPeriodic(GrapheneFluid2D& graphene){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();	
	
	for (int i=0; i<Nx; i++){
		int Bottom=i; //i+0*Nx
		int Top=i+(Ny-1)*Nx;
		graphene.den[Bottom] = graphene.den[Top-Nx];
		graphene.flxX[Bottom] = graphene.flxX[Top-Nx];
		graphene.flxY[Bottom] = graphene.flxY[Top-Nx];
		graphene.den[Top] = graphene.den[Bottom+Nx];
		graphene.flxX[Top] = graphene.flxX[Bottom+Nx];
		graphene.flxY[Top] = graphene.flxY[Bottom+Nx];
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

void BoundaryCondition::Dirichlet::MassFluxX(GrapheneFluid2D& graphene, float L, float R, float T, float B){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int j=0; j<Ny; j++){
		graphene.flxX[0+j*Nx] = L;
		graphene.flxX[Nx-1+j*Nx] = R;
	}
	for (int i=0; i<Nx; i++){
		graphene.flxX[i+(Ny-1)*Nx] = T;
		graphene.flxX[i+0*Nx] = B;
	}
} 
void BoundaryCondition::Dirichlet::MassFluxY(GrapheneFluid2D& graphene, float L, float R, float T, float B){
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	for (int j=0; j<Ny; j++){
		graphene.flxY[0+j*Nx] = L;
		graphene.flxY[Nx-1+j*Nx] = R;
	}
	for (int i=0; i<Nx; i++){
		graphene.flxY[i+(Ny-1)*Nx] = T;
		graphene.flxY[i+0*Nx] = B;
	}
} 

void BoundaryCondition::DyakonovShur::X(GrapheneFluid1D& graphene) {
	int Nx=graphene.SizeX();
	graphene.den[0] = 1.0f;
	graphene.vel[0] = graphene.vel[1];
	graphene.den[Nx-1] = graphene.den[Nx-2];
	graphene.vel[Nx-1] = 1.0f/graphene.den[Nx-1];
}

void BoundaryCondition::DyakonovShur::X(GrapheneFluid2D& graphene) {
	int Nx=graphene.SizeX();
	int Ny=graphene.SizeY();
	 
	for(int j=0;j<Ny;j++){	
		graphene.den[0+j*Nx]=1.0f;			//constant density at x=0
		graphene.den[Nx-1+j*Nx]=graphene.den[Nx-2+j*Nx]; 			//free density at x=L
		graphene.flxX[0+j*Nx] = graphene.flxX[1+j*Nx]*pow(graphene.den[1+j*Nx],-1.5f);			//free flux at x=0
		graphene.flxX[Nx-1+j*Nx] = sqrt(graphene.den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
		graphene.flxY[0+j*Nx] = 0.0f; 					//flux only on x at x=0
		graphene.flxY[Nx-1+j*Nx] = 0.0f ;					//idem at x=L
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
