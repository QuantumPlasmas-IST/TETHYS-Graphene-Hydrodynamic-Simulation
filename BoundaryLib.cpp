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
	int nx=graphene.SizeX();
	graphene.Den[0] = graphene.Den[1];
	graphene.Den[nx - 1] =  graphene.Den[nx - 2];
	graphene.Vel[0] = graphene.Vel[1];
	graphene.Vel[nx - 1] = graphene.Vel[nx - 2];
}
void BoundaryCondition::XPeriodic(GrapheneFluid1D& graphene){
	int nx=graphene.SizeX();
	graphene.Den[0] = graphene.Den[nx - 2];
	graphene.Den[nx - 1] = graphene.Den[1];
	graphene.Vel[0] = graphene.Vel[nx - 2];
	graphene.Vel[nx - 1] = graphene.Vel[1];
}
void BoundaryCondition::XFree(GrapheneFluid2D& graphene){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		graphene.Den[left]=graphene.Den[left + 1];
		graphene.Den[right]=graphene.Den[right - 1];			//free density at x=L
		graphene.FlxY[left] = 0.0f; 					//flux only on x at x=0
		graphene.FlxY[right] = 0.0f ;					//idem at x=L
		graphene.FlxX[left] = graphene.FlxX[left + 1] * pow(graphene.Den[left + 1], -1.5f);			//free flux at x=0
		graphene.FlxX[right] =  graphene.FlxX[right - 1];
	}	
}
void BoundaryCondition::XPeriodic(GrapheneFluid2D& graphene){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		graphene.Den[left]=graphene.Den[right - 1];
		graphene.Den[right]=graphene.Den[1 + j * nx];
		graphene.FlxY[left] = 0.0; 					//flux only on x at x=0
		graphene.FlxY[right] = 0.0 ;					//idem at x=L
		graphene.FlxX[left] = graphene.FlxX[right - 1];
		graphene.FlxX[right] =  graphene.FlxX[left + 1];
	}	
}
void BoundaryCondition::YFree(GrapheneFluid2D& graphene){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		graphene.Den[bottom] = graphene.Den[bottom + nx];
		graphene.FlxX[bottom] = graphene.FlxX[bottom + nx];
		graphene.FlxY[bottom] = graphene.FlxY[bottom + nx];
		graphene.Den[top] = graphene.Den[top - nx];
		graphene.FlxX[top] = graphene.FlxX[top - nx];
		graphene.FlxY[top] = graphene.FlxY[top - nx];
	}	 	
}
void BoundaryCondition::YPeriodic(GrapheneFluid2D& graphene){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		graphene.Den[bottom] = graphene.Den[top - nx];
		graphene.FlxX[bottom] = graphene.FlxX[top - nx];
		graphene.FlxY[bottom] = graphene.FlxY[top - nx];
		graphene.Den[top] = graphene.Den[bottom + nx];
		graphene.FlxX[top] = graphene.FlxX[bottom + nx];
		graphene.FlxY[top] = graphene.FlxY[bottom + nx];
	}
}



void BoundaryCondition::Dirichlet::Density(GrapheneFluid1D& graphene, float left, float right){
	int nx=graphene.SizeX();
	graphene.Den[0] = left;
	graphene.Den[nx - 1] = right;
}
void BoundaryCondition::Dirichlet::VelocityX(GrapheneFluid1D& graphene, float left, float right){
	int nx=graphene.SizeX();
	graphene.Vel[0] = left;
	graphene.Vel[nx - 1] = right;
}  
void BoundaryCondition::Dirichlet::Density(GrapheneFluid2D& graphene, float left, float right, float top, float bottom){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	for (int j=0; j < ny; j++){
		graphene.Den[0 + j * nx] = left;
		graphene.Den[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		graphene.Den[i + (ny - 1) * nx] = top;
		graphene.Den[i + 0 * nx] = bottom;
	}
} 

void BoundaryCondition::Dirichlet::MassFluxX(GrapheneFluid2D& graphene, float left, float right, float top, float bottom){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	for (int j=0; j < ny; j++){
		graphene.FlxX[0 + j * nx] = left;
		graphene.FlxX[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		graphene.FlxX[i + (ny - 1) * nx] = top;
		graphene.FlxX[i + 0 * nx] = bottom;
	}
} 
void BoundaryCondition::Dirichlet::MassFluxY(GrapheneFluid2D& graphene, float left, float right, float top, float bottom){
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	for (int j=0; j < ny; j++){
		graphene.FlxY[0 + j * nx] = left;
		graphene.FlxY[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		graphene.FlxY[i + (ny - 1) * nx] = top;
		graphene.FlxY[i + 0 * nx] = bottom;
	}
} 

void BoundaryCondition::DyakonovShur::X(GrapheneFluid1D& graphene) {
	int nx=graphene.SizeX();
	graphene.Den[0] = 1.0f;
	graphene.Vel[0] = graphene.Vel[1];
	graphene.Den[nx - 1] = graphene.Den[nx - 2];
	graphene.Vel[nx - 1] = 1.0f / graphene.Den[nx - 1];
}

void BoundaryCondition::DyakonovShur::X(GrapheneFluid2D& graphene) {
	int nx=graphene.SizeX();
	int ny=graphene.SizeY();
	 
	for(int j=0; j < ny; j++){
		graphene.Den[0 + j * nx]=1.0f;			//constant density at x=0
		graphene.Den[nx - 1 + j * nx]=graphene.Den[nx - 2 + j * nx]; 			//free density at x=L
		graphene.FlxX[0 + j * nx] = graphene.FlxX[1 + j * nx] * pow(graphene.Den[1 + j * nx], -1.5f);			//free flux at x=0
		graphene.FlxX[nx - 1 + j * nx] = sqrt(graphene.Den[nx - 1 + j * nx]);	//constant current at x=L (flux equals mass)
		graphene.FlxY[0 + j * nx] = 0.0f; 					//flux only on x at x=0
		graphene.FlxY[nx - 1 + j * nx] = 0.0f ;					//idem at x=L
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
