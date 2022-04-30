/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/BoundaryLib.h"
#include "includes/DirichletBoundaryLib.h"

void DirichletBoundaryCondition::Density(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = left;
	fluid_class.Den[nx - 1] = right;

	fluid_class.Umain[0].n()=left;
	fluid_class.Umain[nx-1].n()=right;
	fluid_class.Uaux[0].n()=left;
	fluid_class.Uaux[nx-1].n()=right;
}

void DirichletBoundaryCondition::DensityLeft(Fluid1D& fluid_class, float left){
	int nx=fluid_class.SizeX();
	fluid_class.Den[0] = left;

	fluid_class.Umain[0].n()=left;
	fluid_class.Uaux[0].n()=left;
}

void DirichletBoundaryCondition::DensityRight(Fluid1D& fluid_class, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Den[nx - 1] = right;

	fluid_class.Umain[nx-1].n()=right;
	fluid_class.Uaux[nx-1].n()=right;
}

void DirichletBoundaryCondition::VelocityX(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Vel[0] = left;
	fluid_class.Vel[nx - 1] = right;

	fluid_class.Umain[0].v()=left;
	fluid_class.Umain[nx-1].v()=right;
	fluid_class.Uaux[0].v()=left;
	fluid_class.Uaux[nx-1].v()=right;
}

void DirichletBoundaryCondition::VelocityXLeft(Fluid1D& fluid_class, float left){
	int nx=fluid_class.SizeX();
	fluid_class.Vel[0] = left;

	fluid_class.Umain[0].v()=left;
	fluid_class.Uaux[0].v()=left;
}

void DirichletBoundaryCondition::VelocityXRight(Fluid1D& fluid_class, float right){
	int nx=fluid_class.SizeX();
	fluid_class.Vel[nx - 1] = right;

	fluid_class.Umain[nx-1].v()=right;
	fluid_class.Uaux[nx-1].v()=right;
}

void DirichletBoundaryCondition::Density(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[j * nx] = left;
		fluid_class.Den[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.Den[i + (ny - 1) * nx] = top;
		fluid_class.Den[i] = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[j * nx] = left;
		fluid_class.FlxX[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i + (ny - 1) * nx] = top;
		fluid_class.FlxX[i] = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[j * nx] = left;
		fluid_class.FlxY[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i + (ny - 1) * nx] = top;
		fluid_class.FlxY[i ] = bottom;
	}
}

void DirichletBoundaryCondition::Jet(Fluid2D &fluid_class, float left, float left_width, float right, float right_width) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	int n_width_left= static_cast<int>(static_cast<float>(ny) * left_width);
	int n_width_right= static_cast<int>(static_cast<float>(ny) * right_width);
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

void DirichletBoundaryCondition::DensityRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::MassFluxXRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::MassFluxYRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::DensityLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Den[j * nx] = left;
	}
}

void DirichletBoundaryCondition::MassFluxXLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxX[j * nx] = left;
	}
}

void DirichletBoundaryCondition::MassFluxYLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.FlxY[j * nx] = left;
	}
}

void DirichletBoundaryCondition::DensityTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.Den[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::MassFluxXTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::MassFluxYTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::DensityBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.Den[i] = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxXBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.FlxX[i] = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxYBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.FlxY[i] = bottom;
	}
}

void DirichletBoundaryCondition::Temperature(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Tmp[j * nx] = left;
		fluid_class.Tmp[nx - 1 + j * nx] = right;
	}
	for (int i=0; i < nx; i++){
		fluid_class.Tmp[i + (ny - 1) * nx] = top;
		fluid_class.Tmp[i] = bottom;
	}
}

void DirichletBoundaryCondition::TemperatureRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Tmp[nx - 1 + j * nx] = right;
	}
}

void DirichletBoundaryCondition::TemperatureLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		fluid_class.Tmp[j * nx] = left;
	}
}

void DirichletBoundaryCondition::TemperatureTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		fluid_class.Tmp[i + (ny - 1) * nx] = top;
	}
}

void DirichletBoundaryCondition::TemperatureBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		fluid_class.Tmp[i] = bottom;
	}
}



