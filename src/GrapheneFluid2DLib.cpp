/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
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


float GrapheneFluid2D::DensityFluxX(StateVec2D U) {
	//return U.px() / sqrt(U.n());
	float den=U.n();
	float px=U.px();
	return px / sqrt(den);
}


float GrapheneFluid2D::DensityFluxY(StateVec2D U) {
//return U.py()/ sqrt(U.n());
	float den=U.n();
	float py=U.py();
	return py / sqrt(den);
}

float GrapheneFluid2D::XMomentumFluxX(StateVec2D U) {
	//float 	mass=DensityToMass(U.n());
	//return U.px()*U.px()/mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * U.S()* U.S()* U.n()* U.n();

	float sound=vel_snd;
	float den=U.n();
	float px=U.px();
	float mass=DensityToMass(den);
	float Vxy=U.dxvy();

	return px * px / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den - odd_vis*Vxy;
}

float GrapheneFluid2D::XMomentumFluxY(StateVec2D U) {
	float den=U.n();
	float px=U.px();
	float py=U.py();
	float mass=DensityToMass(den);
	float Vyy=U.dyvy();

	return px * py / mass - odd_vis*Vyy;
}
float GrapheneFluid2D::YMomentumFluxY(StateVec2D U) {
	//float 	mass=DensityToMass(U.n());
	//return U.py()*U.py()/mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * U.S()*U.S()*U.n()*U.n();

	float sound=vel_snd;
	float den=U.n();
	float py=U.py();
	float mass=DensityToMass(den);
	float Vyx=U.dyvx();

	return py * py / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den + odd_vis*Vyx;
}

float GrapheneFluid2D::YMomentumFluxX(StateVec2D U) {
	//float sound=U.S();
	float den=U.n();
	float px=U.px();
	float py=U.py();
	float Vxx=U.dxvx();
	float mass=DensityToMass(den);
	return px * py / mass + odd_vis*Vxx;
}


float GrapheneFluid2D::DensitySource(StateVec2D U) {
	return 0;
}

float GrapheneFluid2D::XMomentumSource(StateVec2D U) {
	return -1.0f*col_freq*U.px()  - cyc_freq*U.py()/sqrt(U.n());
}

float GrapheneFluid2D::YMomentumSource(StateVec2D U) {
	return -1.0f*col_freq*U.py()  + cyc_freq*U.px()/sqrt(U.n());
}

float GrapheneFluid2D::TemperatureSource(StateVec2D U) {
	return 0;
}


GrapheneFluid2D::~GrapheneFluid2D(){
delete[] Umain;
delete[] Umid;
delete[] Den;
delete[] VelX;
delete[] VelY;
delete[] CurX;
delete[] CurY;
delete[] Tmp;
}

float GrapheneFluid2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}



