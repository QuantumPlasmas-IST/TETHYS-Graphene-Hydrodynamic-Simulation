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

void GrapheneFluid2D::MassFluxToVelocity(string grid) {
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



float GrapheneFluid2D::DensityFluxX(GridPoint p, char side) {
	float * den_ptr;
	float * px_ptr;
	float den = 1.0f;
	float px = 0.0f;
	if(p.IsMidGrid){
		den_ptr = Den; // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		px_ptr = FlxX;
	}else{
		den_ptr = den_mid;
		px_ptr = flxX_mid;
	}
	if (side == 'E'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.SE]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
	}
	if (side == 'W'){
		den = 0.5f*(den_ptr[p.NW] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
	}
	return px / sqrt(den);
}

float GrapheneFluid2D::DensityFluxY(GridPoint p, char side) {
	float * den_ptr;
	float * py_ptr;
	float den = 1.0f;
	float py = 0.0f;
	if(p.IsMidGrid){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		den_ptr = Den;
		py_ptr = FlxY;
	}else{
		den_ptr = den_mid;
		py_ptr = flxY_mid;
	}
	if (side == 'N'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.NW]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
	}
	if (side == 'S'){
		den = 0.5f*(den_ptr[p.SE] + den_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
	}
	return py / sqrt(den);
}

float GrapheneFluid2D::XMomentumFluxX(GridPoint p, char side) {
	float * vel_ptr;
	float * den_ptr;
	float * px_ptr;
	float *dvel_ptr;
	float sound =0.0f;
	float den =1.0f;
	float px =0.0f;
	float dvy=0.0f;
	float mass;
	if(p.IsMidGrid){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		vel_ptr = vel_snd_arr;
		den_ptr = Den;
		px_ptr = FlxX;
		dvel_ptr = velY_dx;
	}else{
		vel_ptr = vel_snd_arr_mid;
		den_ptr = den_mid;
		px_ptr = flxX_mid;
		dvel_ptr = velY_dx_mid;
	}
	if (side == 'E'){
		sound= 0.5f*(vel_ptr[p.NE] + vel_ptr[p.SE]);
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.SE]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
		dvy =  0.5f*(dvel_ptr[p.NE] + dvel_ptr[p.SE]);
	}
	if (side == 'W'){
		sound = 0.5f*(vel_ptr[p.NW] + vel_ptr[p.SW]);
		den = 0.5f*(den_ptr[p.NW] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
		dvy =  0.5f*(dvel_ptr[p.NW] + dvel_ptr[p.SW]);
	}
	mass=DensityToMass(den);
	return px * px / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den - odd_vis*dvy;
}

float GrapheneFluid2D::XMomentumFluxY(GridPoint p, char side) {
	float * den_ptr;
	float * px_ptr;
	float * py_ptr;
	float *dvel_ptr;
	float den =1.0f;
	float px =0.0f;
	float py =0.0f;
	float dvy=0.0f;
	float mass;
	if(p.IsMidGrid){
		den_ptr = Den;
		px_ptr = FlxX;
		py_ptr = FlxY;
		dvel_ptr = velY_dy;
	}else{
		den_ptr = den_mid;
		px_ptr = flxX_mid;
		py_ptr = flxY_mid;
		dvel_ptr = velY_dy_mid;
	}
	if (side == 'N'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.NW]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.NW]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
		dvy =  0.5f*(dvel_ptr[p.NE] + dvel_ptr[p.NW]);
	}
	if (side == 'S'){
		den = 0.5f*(den_ptr[p.SE] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.SE] + px_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
		dvy =  0.5f*(dvel_ptr[p.SE] + dvel_ptr[p.SW]);
	}
	mass=DensityToMass(den);
	return px * py / mass - odd_vis*dvy;
}


float GrapheneFluid2D::YMomentumFluxY(GridPoint p, char side) {
	float * vel_ptr;
	float * den_ptr;
	float * py_ptr;
	float * dvel_ptr;
	float sound =0.0f;
	float den =1.0f;
	float py =0.0f;
	float dvx =0.0f;
	float mass;
	if(p.IsMidGrid){
		vel_ptr = vel_snd_arr;
		den_ptr = Den;
		py_ptr = FlxY;
		dvel_ptr = velX_dy;
	}else{
		vel_ptr = vel_snd_arr_mid;
		den_ptr = den_mid;
		py_ptr = flxY_mid;
		dvel_ptr = velX_dy_mid;
	}
	if (side == 'N'){
		sound= 0.5f*(vel_ptr[p.NE] + vel_ptr[p.NW]);
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.NW]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
		dvx =  0.5f*(dvel_ptr[p.NE] + dvel_ptr[p.NW]);
	}
	if (side == 'S'){
		sound = 0.5f*(vel_ptr[p.SE] + vel_ptr[p.SW]);
		den = 0.5f*(den_ptr[p.SE] + den_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
		dvx =  0.5f*(dvel_ptr[p.SE] + dvel_ptr[p.SW]);
	}
	mass=DensityToMass(den);
	return py * py / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den + odd_vis*dvx;
}

float GrapheneFluid2D::YMomentumFluxX(GridPoint p, char side) {
	float * den_ptr;
	float * px_ptr;
	float * py_ptr;
	float * dvel_ptr;
	float den =1.0f;
	float px =0.0f;
	float py =0.0f;
	float dvx=0.0f;
	float mass;
	if(p.IsMidGrid){
		den_ptr = Den;
		px_ptr = FlxX;
		py_ptr = FlxY;
		dvel_ptr = velX_dx;
	}else{
		den_ptr = den_mid;
		px_ptr = flxX_mid;
		py_ptr = flxY_mid;
		dvel_ptr = velX_dx_mid;
	}
	if (side == 'E'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.SE]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.SE]);
		dvx =  0.5f*(dvel_ptr[p.NE] + dvel_ptr[p.NW]);
	}
	if (side == 'W'){
		den = 0.5f*(den_ptr[p.NW] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.NW] + py_ptr[p.SW]);
		dvx =  0.5f*(dvel_ptr[p.NW] + dvel_ptr[p.SW]);
	}
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
delete[] flxX_mid;
delete[] flxY_mid;
delete[] lap_flxX;
delete[] lap_flxY;
delete[] vel_snd_arr;
delete[] vel_snd_arr_mid;
delete[] velX_dx;
delete[] velX_dx_mid;
delete[] velX_dy;
delete[] velX_dy_mid;
delete[] velY_dx;
delete[] velY_dx_mid;
delete[] velY_dy;
delete[] velY_dy_mid;
}

float GrapheneFluid2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}

float
GrapheneFluid2D::TemperatureSource(float n, float flx_x, float flx_y, float den_grad_x, float den_grad_y, float mass, float s) {
	return vel_snd * vel_snd * (den_grad_x * flx_x / sqrt(n) + den_grad_y * flx_y / sqrt(n)  ) / (vel_fer * vel_fer);
}


