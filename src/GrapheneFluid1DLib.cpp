/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/GrapheneFluid1DLib.h"



GrapheneFluid1D::GrapheneFluid1D(SetUpParameters &input_parameters) : Fluid1D(input_parameters) {
	vel_fer = input_parameters.FermiVelocity;
	col_freq = input_parameters.CollisionFrequency;

	param = {vel_snd,vel_fer,0.0f,kin_vis,0.0f,therm_diff,col_freq,0.0f};

	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}

GrapheneFluid1D::~GrapheneFluid1D(){
	delete Den;
	delete Vel ;
	delete Cur ;
	delete Umain;
	delete Umid;
	delete Uaux;
	delete vel_snd_arr ;
	delete GradVel ;
}



float GrapheneFluid1D::VelocityFlux(StateVec U) {

	return 0.25f * U.v() * U.v() + vel_fer * vel_fer * 0.5f * log(U.n()+1.0E-6) + 2.0f * U.S() * U.S() * sqrt(U.n()); //TODO falta o termo dv para a voscosidade
}


float GrapheneFluid1D::DensityFlux(StateVec U) {
	return U.n()*U.v();
}

float GrapheneFluid1D::DensitySource(StateVec U){
	return 0.0f;
}

float GrapheneFluid1D::VelocitySource(StateVec U) {
	return -1.0f * col_freq * U.v();
}






float GrapheneFluid1D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0.0f;
}

float GrapheneFluid1D::VelocitySource(float n, float v, float s, float d3den) {
	return -1.0f * col_freq * v;
}

void GrapheneFluid1D::CflCondition(){
	dx = lengX / ( float ) ( Nx - 1 );
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
}

float GrapheneFluid1D::JacobianSpectralRadius(StateVec U) {
	float SQRT = sqrt(16.0f*sqrt(U.n())*vel_snd*vel_snd + U.v()*U.v() + 8.0f*vel_fer*vel_fer);
	float l1 = abs(3.0f*U.v() + SQRT);
	float l2 = abs(3.0f*U.v() - SQRT);
	return 0.25f*max(l1,l2);
}
