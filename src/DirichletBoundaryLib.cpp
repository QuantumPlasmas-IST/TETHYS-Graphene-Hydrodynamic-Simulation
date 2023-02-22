/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/BoundaryLib.h"
#include "includes/DirichletBoundaryLib.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// ONE DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void DirichletBoundaryCondition::Density(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
//	fluid_class.Den[0] = left;
//	fluid_class.Den[nx - 1] = right;

	fluid_class.Umain[0].n()=left;
	fluid_class.Umain[nx-1].n()=right;
	fluid_class.Uaux[0].n()=left;
	fluid_class.Uaux[nx-1].n()=right;
}

void DirichletBoundaryCondition::DensityLeft(Fluid1D& fluid_class, float left){
//	int nx=fluid_class.SizeX();
//	fluid_class.Den[0] = left;

	fluid_class.Umain[0].n()=left;
	fluid_class.Uaux[0].n()=left;
}

void DirichletBoundaryCondition::DensityRight(Fluid1D& fluid_class, float right){
	int nx=fluid_class.SizeX();
//	fluid_class.Den[nx - 1] = right;

	fluid_class.Umain[nx-1].n()=right;
	fluid_class.Uaux[nx-1].n()=right;
}

void DirichletBoundaryCondition::VelocityX(Fluid1D& fluid_class, float left, float right){
	int nx=fluid_class.SizeX();
//	fluid_class.Vel[0] = left;
//	fluid_class.Vel[nx - 1] = right;

	fluid_class.Umain[0].v()=left;
	fluid_class.Umain[nx-1].v()=right;
	fluid_class.Uaux[0].v()=left;
	fluid_class.Uaux[nx-1].v()=right;
}

void DirichletBoundaryCondition::VelocityXLeft(Fluid1D& fluid_class, float left){
//	int nx=fluid_class.SizeX();
//	fluid_class.Vel[0] = left;

	fluid_class.Umain[0].v()=left;
	fluid_class.Uaux[0].v()=left;
}

void DirichletBoundaryCondition::VelocityXRight(Fluid1D& fluid_class, float right){
	int nx=fluid_class.SizeX();
//	fluid_class.Vel[nx - 1] = right;

	fluid_class.Umain[nx-1].v()=right;
	fluid_class.Uaux[nx-1].v()=right;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// TWO DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void DirichletBoundaryCondition::Density(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Den[j * nx] = left;
		//fluid_class.Den[nx - 1 + j * nx] = right;
		fluid_class.Umain[j * nx].n() = left;
		fluid_class.Umain[nx - 1 + j * nx].n() = right;
	}
	for (int i=0; i < nx; i++){
		//fluid_class.Den[i + (ny - 1) * nx] = top;
		//fluid_class.Den[i] = bottom;
		fluid_class.Umain[i + (ny - 1) * nx].n() = top;
		fluid_class.Umain[i].n() = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxX[j * nx] = left;
		//fluid_class.FlxX[nx - 1 + j * nx] = right;

		fluid_class.Umain[j * nx].px() = left;
		fluid_class.Umain[nx - 1 + j * nx].px() = right;
	}
	for (int i=0; i < nx; i++){
		//fluid_class.FlxX[i + (ny - 1) * nx] = top;
		//fluid_class.FlxX[i] = bottom;

		fluid_class.Umain[i + (ny - 1) * nx].px() = top;
		fluid_class.Umain[i].px() = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxY[j * nx] = left;
		//fluid_class.FlxY[nx - 1 + j * nx] = right;
		fluid_class.Umain[j * nx].py() = left;
		fluid_class.Umain[nx - 1 + j * nx].py() = right;
	}
	for (int i=0; i < nx; i++){
		//fluid_class.FlxY[i + (ny - 1) * nx] = top;
		//fluid_class.FlxY[i ] = bottom;
		fluid_class.Umain[i + (ny - 1) * nx].py() = top;
		fluid_class.Umain[i].py() = bottom;
	}
}

void DirichletBoundaryCondition::Jet(Fluid2D &fluid_class, float left, float left_width, float right, float right_width) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	int n_width_left= static_cast<int>(static_cast<float>(ny) * left_width);
	int n_width_right= static_cast<int>(static_cast<float>(ny) * right_width);
	for (int j=0; j < ny; j++){
		if( j>=(ny-n_width_left)/2 && j<= (ny+n_width_left)/2){
			//fluid_class.FlxX[0 + j * nx] = left;
			fluid_class.Umain[0 + j * nx].px() = left;
		} else{
			//fluid_class.FlxX[0 + j * nx] = 0.0f;
			fluid_class.Umain[0 + j * nx].px() =0.0f;
		}
		if( j>=(ny-n_width_right)/2 && j<= (ny+n_width_right)/2){
			//fluid_class.FlxX[nx - 1 + j * nx] = right;
			fluid_class.Umain[nx - 1 + j * nx].px() = right;
		} else{
			//fluid_class.FlxX[nx - 1 + j * nx] = 0.0f;
			fluid_class.Umain[nx - 1 + j * nx].px() =0.0f;
		}
	}
}

void DirichletBoundaryCondition::DensityRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Den[nx - 1 + j * nx] = right;
		fluid_class.Umain[nx - 1 + j * nx].n() = right;
	}
}

void DirichletBoundaryCondition::MassFluxXRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxX[nx - 1 + j * nx] = right;
		fluid_class.Umain[nx - 1 + j * nx].px() = right;
	}
}

void DirichletBoundaryCondition::MassFluxYRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxY[nx - 1 + j * nx] = right;
		fluid_class.Umain[nx - 1 + j * nx].py() = right;
	}
}

void DirichletBoundaryCondition::DensityLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Den[j * nx] = left;
		fluid_class.Umain[j * nx].n() = left;
	}
}

void DirichletBoundaryCondition::MassFluxXLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxX[j * nx] = left;
		fluid_class.Umain[j * nx].px() = left;
	}
}

void DirichletBoundaryCondition::MassFluxYLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.FlxY[j * nx] = left;
		fluid_class.Umain[j * nx].py() = left;
	}
}

void DirichletBoundaryCondition::DensityTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		//fluid_class.Den[i + (ny - 1) * nx] = top;
		fluid_class.Umain[i + (ny - 1) * nx].n() = top;
	}
}

void DirichletBoundaryCondition::MassFluxXTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		//fluid_class.FlxX[i + (ny - 1) * nx] = top;
		fluid_class.Umain[i + (ny - 1) * nx].px() = top;
	}
}

void DirichletBoundaryCondition::MassFluxYTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		//fluid_class.FlxY[i + (ny - 1) * nx] = top;
		fluid_class.Umain[i + (ny - 1) * nx].py() = top;
	}
}

void DirichletBoundaryCondition::DensityBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		//fluid_class.Den[i] = bottom;
		fluid_class.Umain[i].n() = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxXBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		//fluid_class.FlxX[i] = bottom;
		fluid_class.Umain[i].px() = bottom;
	}
}

void DirichletBoundaryCondition::MassFluxYBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		//fluid_class.FlxY[i] = bottom;
		fluid_class.Umain[i].py() = bottom;
	}
}

void DirichletBoundaryCondition::Temperature(Fluid2D& fluid_class, float left, float right, float top, float bottom){
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Tmp[j * nx] = left;
		//fluid_class.Tmp[nx - 1 + j * nx] = right;
		fluid_class.Umain[j * nx].tmp() = left;
		fluid_class.Umain[nx - 1 + j * nx].tmp() = right;

	}
	for (int i=0; i < nx; i++){
		//fluid_class.Tmp[i + (ny - 1) * nx] = top;
		//fluid_class.Tmp[i] = bottom;
		fluid_class.Umain[i + (ny - 1) * nx].tmp() = top;
		fluid_class.Umain[i].tmp() = bottom;
	}
}

void DirichletBoundaryCondition::TemperatureRight(Fluid2D &fluid_class, float right) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Tmp[nx - 1 + j * nx] = right;
		fluid_class.Umain[nx - 1 + j * nx].tmp() = right;
	}
}

void DirichletBoundaryCondition::TemperatureLeft(Fluid2D &fluid_class, float left) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int j=0; j < ny; j++){
		//fluid_class.Tmp[j * nx] = left;
		fluid_class.Umain[j * nx].tmp() = left;
	}
}

void DirichletBoundaryCondition::TemperatureTop(Fluid2D &fluid_class, float top) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	for (int i=0; i < nx; i++){
		//fluid_class.Tmp[i + (ny - 1) * nx] = top;
		fluid_class.Umain[i + (ny - 1) * nx].tmp() = top;
	}
}

void DirichletBoundaryCondition::TemperatureBottom(Fluid2D &fluid_class, float bottom) {
	int nx=fluid_class.SizeX();
	for (int i=0; i < nx; i++){
		//fluid_class.Tmp[i] = bottom;
		fluid_class.Umain[i].tmp() = bottom;
	}
}



