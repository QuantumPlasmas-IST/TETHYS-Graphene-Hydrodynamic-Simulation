/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
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



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// ONE-DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



void BoundaryCondition::XFree(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Umain[0] =fluid_class.Umain[1];
	fluid_class.Uaux[0] =fluid_class.Uaux[1];
	fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];
}

void BoundaryCondition::XFree(DiracGraphene1D& fluid_class){
    int nx=fluid_class.SizeX();
    fluid_class.Umain[0] =fluid_class.Umain[1];
    fluid_class.Uaux[0] =fluid_class.Uaux[1];
    fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
    fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];

    fluid_class.HoleUmain[0] =fluid_class.HoleUmain[1];
    fluid_class.HoleUmain[nx-1] =fluid_class.HoleUmain[nx-2];
}

void BoundaryCondition::XFreeLeft(Fluid1D &fluid_class) {
	fluid_class.Umain[0] =fluid_class.Umain[1];
	fluid_class.Uaux[0] =fluid_class.Uaux[1];
}

void BoundaryCondition::XFreeLeft(DiracGraphene1D &fluid_class) {
    fluid_class.Umain[0] =fluid_class.Umain[1];
    fluid_class.Uaux[0] =fluid_class.Uaux[1];

    fluid_class.HoleUmain[0] =fluid_class.HoleUmain[1];
}

void BoundaryCondition::XFreeRight(Fluid1D &fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];
}

void BoundaryCondition::XFreeRight(DiracGraphene1D &fluid_class) {
    int nx=fluid_class.SizeX();
    fluid_class.Umain[nx-1] =fluid_class.Umain[nx-2];
    fluid_class.Uaux[nx-1] =fluid_class.Uaux[nx-2];

    fluid_class.HoleUmain[nx-1] =fluid_class.HoleUmain[nx-2];
}

void BoundaryCondition::XPeriodic(Fluid1D& fluid_class){
	int nx=fluid_class.SizeX();
	fluid_class.Umain[0] =fluid_class.Umain[nx-2];
	fluid_class.Uaux[0] =fluid_class.Uaux[nx-2];
	fluid_class.Umain[nx-1] =fluid_class.Umain[1];
	fluid_class.Uaux[nx-1] =fluid_class.Uaux[1];
}

void BoundaryCondition::XPeriodic(DiracGraphene1D& fluid_class){
    int nx=fluid_class.SizeX();
    fluid_class.Umain[0] =fluid_class.Umain[nx-2];
    fluid_class.Uaux[0] =fluid_class.Uaux[nx-2];
    fluid_class.Umain[nx-1] =fluid_class.Umain[1];
    fluid_class.Uaux[nx-1] =fluid_class.Uaux[1];

    fluid_class.HoleUmain[0] =fluid_class.HoleUmain[nx-2];
    fluid_class.HoleUmain[nx-1] =fluid_class.HoleUmain[1];
}

//TODO condicoes fronteira para os diracs e ver o que se passa com as periodicas

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// TWO-DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
		fluid_class.Umain[pos] = fluid_class.Umain[neighbour];
	}
}
void BoundaryCondition::XFreeLeft(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int left;
		left= 0 + j * nx;
		fluid_class.Umain[left] = fluid_class.Umain[left +1];
	}
}
void BoundaryCondition::XFreeRight(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int right;
		right = nx - 1 + j * nx;
		fluid_class.Umain[right] = fluid_class.Umain[right-1];
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
		fluid_class.Umain[left] = fluid_class.Umain[right - 1];
		fluid_class.Umain[right] =  fluid_class.Umain[left + 1];

	}	
}


void BoundaryCondition::XPeriodic(DiracGraphene2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for(int j=0; j < ny; j++){
		int left;
		left = 0 + j * nx;
		int right;
		right = nx - 1 + j * nx;
		fluid_class.Umain[left] = fluid_class.Umain[right - 1];
		fluid_class.Umain[right] =  fluid_class.Umain[left + 1];
		fluid_class.HoleUmain[left] = fluid_class.HoleUmain[right - 1];
		fluid_class.HoleUmain[right] =  fluid_class.HoleUmain[left + 1];
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
		fluid_class.Umain[bottom] = fluid_class.Umain[bottom + nx];
		fluid_class.Umain[top] = fluid_class.Umain[top - nx];

	}
}
void BoundaryCondition::YFree(Fluid2D &fluid_class, int y_limit) {
	for (int i=0; i <fluid_class.SizeX(); i++){
		int pos = i + (fluid_class.SizeY() - 1) * fluid_class.SizeX()*y_limit;
		int neighbour=pos+(1-2*y_limit)*fluid_class.SizeX();
		fluid_class.Umain[pos] = fluid_class.Umain[neighbour];
	}
}
void BoundaryCondition::YFreeTop(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Umain[top] = fluid_class.Umain[top - nx];

	}
}
void BoundaryCondition::YFreeBottom(Fluid2D &fluid_class) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		int bottom;
		bottom = i;
		fluid_class.Umain[bottom] = fluid_class.Umain[bottom + nx];
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
		fluid_class.Umain[bottom] = fluid_class.Umain[top - nx];
		fluid_class.Umain[top] = fluid_class.Umain[bottom + nx];
	}
}


void BoundaryCondition::YPeriodic(DiracGraphene2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		int bottom;
		bottom = i; //i+0*nx
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Umain[bottom] = fluid_class.Umain[top - nx];
		fluid_class.Umain[top] = fluid_class.Umain[bottom + nx];
		fluid_class.HoleUmain[bottom] = fluid_class.HoleUmain[top - nx];
		fluid_class.HoleUmain[top] = fluid_class.HoleUmain[bottom + nx];
	}
}

void BoundaryCondition::YClosedFreeSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){

		int bottom;
		bottom = i; //i+0*nx
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Umain[bottom].n() = fluid_class.Umain[bottom+nx].n();
		fluid_class.Umain[top].n() = fluid_class.Umain[top-nx].n();

		fluid_class.Umain[bottom].px() = fluid_class.Umain[top - nx].px();
		fluid_class.Umain[bottom].py() = 0.0f; //Slope * fluid_class.FlxX[top - nx];
		fluid_class.Umain[top].px() = fluid_class.Umain[bottom + nx].px();
		fluid_class.Umain[top].py() = 0.0f; //-1.0f * Slope * fluid_class.FlxX[bottom + nx];

	}
}
void BoundaryCondition::YClosedNoSlip(Fluid2D& fluid_class){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++) {
		int bottom;
		bottom = i; //i+0*nx
		int top;
		top = i + (ny - 1) * nx;
		fluid_class.Umain[bottom].n() = fluid_class.Umain[bottom+nx].n();
		fluid_class.Umain[top].n() = fluid_class.Umain[top-nx].n();
		fluid_class.Umain[bottom].px() = 0.0f;
		fluid_class.Umain[bottom].py() = 0.0f; //Slope * fluid_class.FlxX[top - nx];
		fluid_class.Umain[top].px() = 0.0f;
		fluid_class.Umain[top].py() = 0.0f; //-1.0f * Slope * fluid_class.FlxX[bottom + nx];

	}
}





