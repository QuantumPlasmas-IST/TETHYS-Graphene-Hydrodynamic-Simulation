/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/GrapheneFluid2DLib.h"



GrapheneFluid2D::GrapheneFluid2D(SetUpParameters &input_parameters) : Fluid2D(input_parameters) {
	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency;
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency;
	therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity
	char buffer [100];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
	file_infix = buffer;
}

void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

/*
void GrapheneFluid2D::MassFluxToVelocity(const string grid) {
float den;
	if(grid=="MainGrid"){
		for(int c=0; c <= Nx * Ny - 1; c++){
			den = Den[c];
			VelX[c] = FlxX[c] / sqrt(den*den*den);
			VelY[c] = FlxY[c] / sqrt(den*den*den);
		}
	}
	if(grid=="MidGrid"){
		for(int c=0; c <= (Nx-1) * (Ny-1) - 1; c++){
			den = den_mid[c];
			velX_mid[c] = flxX_mid[c] / sqrt(den*den*den);
			velY_mid[c] = flxY_mid[c] / sqrt(den*den*den);
		}
	}
}
*/

void GrapheneFluid2D::CflCondition(){ // Eventual redefinition
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
	/*  CFL condition for FTCS method
	if(kin_vis>0.0f&& kin_vis*dt > dx*dx*0.25f){
		dt = 0.8f*0.25f*dx*dx/kin_vis;
	}*/
	//  CFL condition for (1,9) Weighted explicit method
	if(kin_vis>0.0f&& kin_vis*dt > dx*dx*0.5f){
		dt = 0.8f*0.5f*dx*dx/kin_vis;
	}
	if(therm_diff>0.0f&& therm_diff*dt > dx*dx*0.5f){
		dt = 0.8f*0.5f*dx*dx/therm_diff;
	}
}



float GrapheneFluid2D::DensityFluxX(GridPoint2D p, char side) {
	float den;
	float px;

	px = SideAverage(ptr_px, p, side);
	den = SideAverage(ptr_den, p, side);

	return px / sqrt(den);
}

float GrapheneFluid2D::DensityFluxY(GridPoint2D p, char side) {
	float den;
	float py;
	py = SideAverage(ptr_py, p, side);
	den = SideAverage(ptr_den, p, side);
	return py / sqrt(den);
}

float GrapheneFluid2D::XMomentumFluxX(GridPoint2D p, char side) {

	float sound;
	float den;
	float px;
	float dvy;
	float mass;
	float d2den;

	sound = SideAverage(ptr_snd, p, side);
	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	dvy = SideAverage(ptr_velYdx, p, side);

	mass=DensityToMass(den);

	d2den = SideAverage(ptr_lap_den, p, side);

	return px * px / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den - odd_vis*dvy   ;
}

float GrapheneFluid2D::XMomentumFluxY(GridPoint2D p, char side) {

	float den;
	float py;
	float px;
	float dvy;
	float mass;



	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	dvy = SideAverage(ptr_velYdy, p, side);
	mass=DensityToMass(den);

	return px * py / mass - odd_vis*dvy;
}


float GrapheneFluid2D::YMomentumFluxY(GridPoint2D p, char side) {

	float sound ;
	float den;
	float py;
	float dvx;
	float mass;

	sound = SideAverage(ptr_snd, p, side);
	den = SideAverage(ptr_den, p, side);
	py = SideAverage(ptr_py, p, side);
	dvx = SideAverage(ptr_velXdy, p, side);
	mass=DensityToMass(den);

	float d2den = SideAverage(ptr_lap_den, p, side);

	return py * py / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den + odd_vis*dvx ;
}

float GrapheneFluid2D::YMomentumFluxX(GridPoint2D p, char side) {

	float den;
	float px;
	float py;
	float dvx;
	float mass;

	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	dvx = SideAverage(ptr_velXdx, p, side);
	mass=DensityToMass(den);

	return px * py / mass  + odd_vis*dvx;
}


float GrapheneFluid2D::DensitySource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

float GrapheneFluid2D::XMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_x  - cyc_freq*flx_y/sqrt(n);
}

float GrapheneFluid2D::YMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_y  + cyc_freq*flx_x/sqrt(n);
}


GrapheneFluid2D::~GrapheneFluid2D(){
delete[] Den;
delete[] VelX;
delete[] VelY;
delete[] FlxX;
delete[] FlxY;
delete[] CurX;
delete[] CurY;
delete[] den_mid;
delete[] velX_mid;
delete[] velY_mid;
delete[] flxX_mid;
delete[] flxY_mid;
delete[] lap_flxX;
delete[] lap_flxY;
delete[] vel_snd_arr;
delete[] vel_snd_arr_mid;
delete[] den_dx;
delete[] den_dy;
delete[] den_dx_mid;
delete[] den_dy_mid;
delete[] velX_dx;
delete[] velX_dx_mid;
delete[] velX_dy;
delete[] velX_dy_mid;
delete[] velY_dx;
delete[] velY_dx_mid;
delete[] velY_dy;
delete[] velY_dy_mid;
delete[] lap_den;
delete[] lap_tmp;
delete[] lap_den_mid;
delete[] tmp_mid;
delete[] Tmp;
}

float GrapheneFluid2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}

float
GrapheneFluid2D::TemperatureSource(float n, float flx_x, float flx_y, float den_grad_x, float den_grad_y, float mass, float s) {
	return vel_snd * vel_snd * (den_grad_x * flx_x / sqrt(n) + den_grad_y * flx_y / sqrt(n)  ) / (vel_fer * vel_fer);
}


