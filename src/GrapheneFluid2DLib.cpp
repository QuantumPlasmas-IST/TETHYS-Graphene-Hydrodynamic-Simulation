//
// Created by pcosme on 27/12/2020.
//

#include "GrapheneFluid2DLib.h"

GrapheneFluid2D::GrapheneFluid2D(SetUpParameters &input_parameters) : Fluid2D(input_parameters) {
	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency;
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency;
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2fwc=%.2f", vel_snd, vel_fer, kin_vis, col_freq,cyc_freq);
	file_infix = buffer;
}

void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

void GrapheneFluid2D::MassFluxToVelocity(){
//#pragma omp parallel for default(none) shared(VelX,VelY,FlxX,FlxY,Den,Nx,Ny)
float den;
	for(int c=0; c <= Nx * Ny - 1; c++){
//		VelX[c]= FlxX[c] * pow(Den[c], -1.5f);
//		VelY[c]= FlxY[c] * pow(Den[c], -1.5f);
den=Den[c];
		VelX[c]= FlxX[c] / sqrt(den*den*den);
		VelY[c]= FlxY[c] / sqrt(den*den*den);
		//CurX[c] = VelX[c] * Den[c];
		//CurY[c] = VelY[c] * Den[c];
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
}

float  GrapheneFluid2D::DensityFluxX(float n, float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused))  float s){
	float f_1;
	f_1 = flx_x / sqrt(n);
	return f_1;
}

float  GrapheneFluid2D::DensityFluxY(float n, __attribute__((unused)) float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused))  float s){
	float f_1;
	f_1 = flx_y / sqrt(n);
	return f_1;
}

float  GrapheneFluid2D::MassFluxXFluxX(float n, float flx_x, __attribute__((unused)) float flx_y, float mass, float s){
	float f_2;
	f_2 = flx_x * flx_x / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_2;
}

float  GrapheneFluid2D::MassFluxXFluxY(__attribute__((unused)) float n, float flx_x, float flx_y, float mass, __attribute__((unused)) float s){
	float f_2;
	f_2 = flx_x * flx_y / mass;
	return f_2;
}

float  GrapheneFluid2D::MassFluxYFluxX(__attribute__((unused)) float n, float flx_x, float flx_y, float mass, __attribute__((unused)) float s){
	float f_3;
	f_3 = flx_x * flx_y / mass;
	return f_3;
}

float  GrapheneFluid2D::MassFluxYFluxY(float n, __attribute__((unused)) float flx_x, float flx_y, float mass, float s){
	float f_3;
	f_3 = flx_y * flx_y / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_3;
}

void GrapheneFluid2D::MagneticSourceSemiAnalytic(){
	float px_0,py_0,sqrtn_0;
	float wc=cyc_freq;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 * cos(wc * dt / sqrtn_0) - py_0 * sin(wc * dt / sqrtn_0);
			FlxY[kp]= px_0 * sin(wc * dt / sqrtn_0) + py_0 * cos(wc * dt / sqrtn_0);
		}
	}
}

/*
void GrapheneFluid2D::MagneticSourceFtcs(){
	float px_0,py_0,sqrtn_0;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 -  dt * cyc_freq * py_0 / sqrtn_0;
			FlxY[kp]= py_0 +  dt * cyc_freq * px_0 / sqrtn_0;
		}
	}
}
*/
float GrapheneFluid2D::DensitySource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

float GrapheneFluid2D::MassFluxXSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_x  - cyc_freq*flx_y/sqrt(n);
}

float GrapheneFluid2D::MassFluxYSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
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
}

