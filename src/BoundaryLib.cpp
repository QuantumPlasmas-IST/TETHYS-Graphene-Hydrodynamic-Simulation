/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/BoundaryLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif


using namespace H5;
using namespace std;


float  BoundaryCondition::Slope=0.0f;

int * BoundaryCondition::BottomEdge=nullptr;
int * BoundaryCondition::TopEdge=nullptr;


void BoundaryCondition::XFree(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[1];
	fluid_class.Den[nx - 1] =  fluid_class.Den[nx - 2];
	fluid_class.Vel[0] = fluid_class.Vel[1];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[nx - 2];

	fluid_class.Umain[0] =fluid_class.Umain[1];
	fluid_class.Uaux[0] =fluid_class.Uaux[1];
	fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];
}

void BoundaryCondition::XFreeLeft(Fluid1D &fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[1];
	fluid_class.Vel[0] = fluid_class.Vel[1];

	fluid_class.Umain[0] =fluid_class.Umain[1];
	fluid_class.Uaux[0] =fluid_class.Uaux[1];
}

void BoundaryCondition::XFreeRight(Fluid1D &fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Den[nx - 1] =  fluid_class.Den[nx - 2];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[nx - 2];

	fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];
}

void BoundaryCondition::XPeriodic(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = fluid_class.Den[nx - 2];
	fluid_class.Den[nx - 1] = fluid_class.Den[1];
	fluid_class.Vel[0] = fluid_class.Vel[nx - 2];
	fluid_class.Vel[nx - 1] = fluid_class.Vel[1];

	fluid_class.Umain[0] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[0] =fluid_class.Uaux[nx-2];
	fluid_class.Umain[nx-1] =fluid_class.Umain[1];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[1];

}
void BoundaryCondition::XFree(Fluid2D& fluid_class){
	XFree(fluid_class, 1);
	XFree(fluid_class, 0);
}
void BoundaryCondition::XFree(Fluid2D &fluid_class, int x_limit) {
	for(int j=0; j < fluid_class.SizeY(); j++){
		int pos;
		int neighbour;
		pos = x_limit * fluid_class.SizeX() + j * fluid_class.SizeX();
		neighbour=pos+(1-2*x_limit);
		fluid_class.Den[pos]=fluid_class.Den[neighbour];
		fluid_class.FlxY[pos] = fluid_class.FlxY[neighbour] ;
		fluid_class.FlxX[pos] = fluid_class.FlxX[neighbour] ;
	}
}
void BoundaryCondition::XFreeLeft(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int left;
		left= 0 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[left + 1];
		fluid_class.FlxY[left] = fluid_class.FlxY[left + 1] ;//flux only on x at x=0
		fluid_class.FlxX[left] = fluid_class.FlxX[left + 1] ;//free flux at x=0
	}
}
void BoundaryCondition::XFreeRight(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int right;
		right = nx - 1 + j * nx;
		fluid_class.Den[right]=fluid_class.Den[right - 1];			//free density at x=L
		fluid_class.FlxY[right] = fluid_class.FlxY[right - 1] ;					//idem at x=L
		fluid_class.FlxX[right] =  fluid_class.FlxX[right - 1];
	}
}

void BoundaryCondition::XPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int left;
		left = 0 + j * nx;
		int right;
		right = nx - 1 + j * nx;
		fluid_class.Den[left]=fluid_class.Den[right - 1];
		fluid_class.Den[right]=fluid_class.Den[1 + j * nx];
		fluid_class.FlxY[left] = 0.0f; 					//flux only on x at x=0
		fluid_class.FlxY[right] = 0.0f ;					//idem at x=L
		fluid_class.FlxX[left] = fluid_class.FlxX[right - 1];
		fluid_class.FlxX[right] =  fluid_class.FlxX[left + 1];
	}	
}
void BoundaryCondition::YFree(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom;
		bottom = i; //i+0*nx
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[top - nx];
	}
}
void BoundaryCondition::YFree(Fluid2D &fluid_class, int y_limit) {
	for (int i=0; i <fluid_class.SizeX(); i++){
		int pos = i + (fluid_class.SizeY() - 1) * fluid_class.SizeX()*y_limit;
		int neighbour=pos+(1-2*y_limit)*fluid_class.SizeX();
		fluid_class.Den[pos]=fluid_class.Den[neighbour];
		fluid_class.FlxY[pos] = fluid_class.FlxY[neighbour] ;
		fluid_class.FlxX[pos] = fluid_class.FlxX[neighbour] ;
	}
}
void BoundaryCondition::YFreeTop(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[top] = fluid_class.FlxY[top - nx];
	}
}
void BoundaryCondition::YFreeBottom(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		int bottom;
		bottom = i; //i+0*nx
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[bottom] = fluid_class.FlxY[bottom + nx];
	}
}
void BoundaryCondition::YPeriodic(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom;
		bottom = i; //i+0*nx
		int top;
		top = i + (ny - 1) * nx;
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
	for (int i=0; i < nx; i++){
		int j_top;
		j_top = TopEdge[i];
		int j_bot;
		j_bot = BottomEdge[i];
		int bottom;
		bottom = i + j_bot * nx;
		int top;
		top = i + j_top * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.Den[top] = fluid_class.Den[top - nx];

		fluid_class.FlxX[bottom] = fluid_class.FlxX[top - nx];
		fluid_class.FlxY[bottom] = Slope * fluid_class.FlxX[top - nx];
		fluid_class.FlxX[top] = fluid_class.FlxX[bottom + nx];
		fluid_class.FlxY[top] = -1.0f * Slope * fluid_class.FlxX[bottom + nx];

	}
}
void BoundaryCondition::YClosedNoSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++) {
		int j_top;
		j_top = TopEdge[i];
		int j_bot;
		j_bot = BottomEdge[i];
		int bottom;
		bottom = i + j_bot * nx; //i+0*nx
		int top;
		top = i + j_top * nx;
		fluid_class.Den[bottom] = fluid_class.Den[bottom + nx];
		fluid_class.FlxX[bottom] = 0.0f;
		fluid_class.FlxY[bottom] = 0.0f;
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxX[top] = 0.0f;
		fluid_class.FlxY[top] = 0.0f;
	}
}

void BoundaryCondition::SetBottomEdge(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	BottomEdge = new int[nx]();
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny/2;j++) {
			if(abs(static_cast<float>(j) - Slope * static_cast<float>(i)) <= 0.5f){
				BottomEdge[i] = j;
			}
		}
	}
}




void BoundaryCondition::SetTopEdge(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	TopEdge = new int[nx]();
	for(int i=0;i<nx;i++){
		for(int j=ny-1;j>ny/2;j--) {
			if(abs(static_cast<float>(j - (ny-1)) + Slope * static_cast<float>(i) ) <= 0.5f){
				TopEdge[i] = j;
			}
		}
	}
}





void BoundaryCondition::SetSlope(float boundary_slope) { Slope=boundary_slope;}

float BoundaryCondition::GetSlope() {return Slope;}


