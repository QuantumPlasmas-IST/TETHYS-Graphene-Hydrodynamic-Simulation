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


void BoundaryCondition::XFree(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[1];
	fluid_class.Den[nx - 1] =  fluid_class.Den[nx - 2];
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[nx - 2];
}
void BoundaryCondition::XPeriodic(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[nx - 2];
	fluid_class.Den[nx - 1] = fluid_class.Den[1];
	fluid_class.Vel[0] = fluid_class.Vel[nx - 2];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[1];
}
void BoundaryCondition::XFree(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[left + 1];
		fluid_class.Den[right]=fluid_class.Den[right - 1];			//free density at x=L
		fluid_class.FlxY[left] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxY[right] = 0.0f ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[left + 1] * pow(fluid_class.Den[left + 1], -1.5f);			//free flux at x=0
		fluid_class.FlxX[right] =  fluid_class.FlxX[right - 1];
	}	
}
void BoundaryCondition::XPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	
	for(int j=0; j < ny; j++){
		int left= 0 + j * nx;
		int right= nx - 1 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[right - 1];
		fluid_class.Den[right]=fluid_class.Den[1 + j * nx];
		fluid_class.FlxY[left] = 0.0; 					//flux only on x at x=0
		fluid_class.FlxY[right] = 0.0 ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[right - 1];
		fluid_class.FlxX[right] =  fluid_class.FlxX[left + 1];
	}	
}
void BoundaryCondition::YFree(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[top - nx];
	}	 	
}
void BoundaryCondition::YPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[top - nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[top - nx];
		fluid_class.Den[top] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[bottom + nx];
	}
}

void BoundaryCondition::YClosedFreeSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[top] =  0.0f;
	}
}
void BoundaryCondition::YClosedNoSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom=i; //i+0*nx
		int top= i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] =  0.0f;
		fluid_class.FlxY[bottom] =  0.0f;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] =  0.0f;
		fluid_class.FlxY[top] =  0.0f;
	}
}

void BoundaryCondition::Dirichlet::Density(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = left;
	fluid_class.Den[nx - 1] = right;
}
void BoundaryCondition::Dirichlet::VelocityX(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Vel[0] = left;
	fluid_class.Vel[nx - 1] = right;
}  
void BoundaryCondition::Dirichlet::Density(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[0 + j * nx] = left;
		fluid_class.Den[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.Den[i + (ny - 1) * nx] = top;
		fluid_class.Den[i + 0 * nx] = bottom;
	}
} 

void BoundaryCondition::Dirichlet::MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[0 + j * nx] = left;
		fluid_class.FlxX[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i + (ny - 1) * nx] = top;
		fluid_class.FlxX[i + 0 * nx] = bottom;
	}
} 
void BoundaryCondition::Dirichlet::MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[0 + j * nx] = left;
		fluid_class.FlxY[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i + (ny - 1) * nx] = top;
		fluid_class.FlxY[i + 0 * nx] = bottom;
	}
}

void BoundaryCondition::Dirichlet::Jet(Fluid2D &fluid_class, float left, float left_width, float right, float right_width) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	int n_width_left= static_cast<int>(ny * left_width);
	int n_width_right= static_cast<int>(ny * right_width);
	for (int j=0; j < ny; j++){
		if( j>=(ny-n_width_left)/2 && j<= (ny+n_width_left)/2){
			fluid_class.FlxX[0 + j * nx] = left;
		} else{
			fluid_class.FlxX[0 + j * nx] = 0.0f;
		}
		if( j>=(ny-n_width_right)/2 && j<= (ny+n_width_right)/2){
			fluid_class.FlxX[nx - 1 + j * nx] = right;
		} else{
			fluid_class.FlxX[nx - 1 + j * nx] = 0.0f;
		}
	}
}

void BoundaryCondition::DyakonovShur::X(GrapheneFluid1D& fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = 1.0f;
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Den[nx - 1] = fluid_class.Den[nx - 2];
	fluid_class.Vel[nx - 1] = 1.0f / fluid_class.Den[nx - 1];
}

void BoundaryCondition::DyakonovShur::X(GrapheneFluid2D& fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	 
	for(int j=0; j < ny; j++){
		fluid_class.Den[0 + j * nx]=1.0f;			//constant density at x=0
		fluid_class.Den[nx - 1 + j * nx]=fluid_class.Den[nx - 2 + j * nx]; 			//free density at x=L
		fluid_class.FlxX[0 + j * nx] = fluid_class.FlxX[1 + j * nx] * pow(fluid_class.Den[1 + j * nx], -1.5f);			//free flux at x=0
		fluid_class.FlxX[nx - 1 + j * nx] = sqrt(fluid_class.Den[nx - 1 + j * nx]);	//constant current at x=L (flux equals mass)
		fluid_class.FlxY[0 + j * nx] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxY[nx - 1 + j * nx] = 0.0f ;					//idem at x=L
	}	
}

void BoundaryCondition::DyakonovShur::YFree(GrapheneFluid2D& fluid_class) {
	BoundaryCondition enclosing;	
	enclosing.YFree(fluid_class);
}

void BoundaryCondition::DyakonovShur::YPeriodic(GrapheneFluid2D& fluid_class) {
	BoundaryCondition enclosing;
	enclosing.YPeriodic(fluid_class);
}


void BoundaryCondition::DyakonovShur::YClosedFreeSlip(GrapheneFluid2D &fluid_class){
	BoundaryCondition enclosing;
	enclosing.YClosedFreeSlip(fluid_class);
}


void BoundaryCondition::DyakonovShur::YClosedNoSlip(GrapheneFluid2D &fluid_class){
	BoundaryCondition enclosing;
	enclosing.YClosedNoSlip(fluid_class);
}