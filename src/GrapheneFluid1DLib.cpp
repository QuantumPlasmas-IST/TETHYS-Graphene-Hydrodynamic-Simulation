/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
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
	delete den_mid ;
	delete vel_mid ;
	delete DenCor ;
	delete VelCor ;
	delete CurCor ;
	delete vel_snd_arr ;
	delete GradVel ;
	delete grad_vel_mid ;
}

/*
float GrapheneFluid1D::DensityFlux(float n,float v,float __attribute__((unused)) s){
	float f_1;
	f_1 = n * v;
	return f_1;
}

float GrapheneFluid1D::VelocityFlux(float n, float v, float dv, float s, float d2n) {
	float f_2;
		f_2 = 0.25f * v * v + vel_fer * vel_fer * 0.5f * log(n) + 2.0f * s * s * sqrt(n)- kin_vis * dv + 0.0f*0.5f*d2n/sqrt(n);
	return f_2;
}
*/
float GrapheneFluid1D::VelocityFlux(GridPoint1D p, char side) {
	float v= SideAverage(ptr_vel,p,side);
	float n=SideAverage(ptr_den,p,side);
	float s= SideAverage(ptr_snd,p,side);
	float dv=SideAverage(ptr_veldx,p,side);
	return 0.25f * v * v + vel_fer * vel_fer * 0.5f * log(n) + 2.0f * s * s * sqrt(n)- kin_vis * dv;
}

float GrapheneFluid1D::VelocityFlux(StateVec1D U) {
	return 0.25f * U.v() * U.v() + vel_fer * vel_fer * 0.5f * log(U.n()) + 2.0f * vel_snd * vel_snd * sqrt(U.n()) ; //TODO falta o termo dv para a voscosidade
	// TODO falta fazer para velocidade do som variavel
}
float GrapheneFluid1D::DensityFlux(GridPoint1D p, char side) {
	float v= SideAverage(ptr_vel,p,side);
	float n=SideAverage(ptr_den,p,side);
	return n * v;
}
float GrapheneFluid1D::DensityFlux(StateVec1D U) {
	return U.n()*U.v();
}








float GrapheneFluid1D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0.0f;
}

float GrapheneFluid1D::VelocitySource(float n, float v, float s, float d3den) {
	return -1.0f * col_freq * v - 0.0f*d3den;
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

float GrapheneFluid1D::JacobianSpectralRadius(StateVec1D U) {
	float SQRT= sqrt(16.0f*sqrt(U.n())*vel_snd*vel_snd+U.v()*U.v()+8.0f*vel_fer*vel_fer   );
	float l1=abs(3.0f*U.v() + SQRT);
	float l2=abs(3.0f*U.v() - SQRT);
	return 0.25f*max(l1,l2);
}




/*
void GrapheneFluid1D::BohmPotencial(string grid) {
	if(grid=="main"){
		for ( int i = 1; i <= Nx-2 ; i++ )
		{
			lap_den[i] = ( Den[i-1]-2.0f*Den[i] + Den[i+1]) / (dx*dx);
		}
		lap_den[0]= (2.0f*Den[0]-5.0f*Den[1]+4.0f*Den[2]-1.0f*Den[3])/ (dx*dx);
		lap_den[Nx-1]=(-2.0f*Den[Nx-1]+5.0f*Den[Nx-2]-4.0f*Den[Nx-3]+1.0f*Den[Nx-4])/ (dx*dx);
	}
	if(grid=="mid"){
		for ( int i = 1; i <= Nx-2 ; i++ )
		{
			lap_den_mid[i] = ( den_mid[i-1]-2.0f*den_mid[i] + den_mid[i+1]) / (dx*dx);
		}
		lap_den_mid[0]= (2.0f*den_mid[0]-5.0f*Den[1]+4.0f*den_mid[2]-1.0f*den_mid[3])/ (dx*dx);
		lap_den_mid[Nx-1]=(-2.0f*den_mid[Nx-1]+5.0f*den_mid[Nx-2]-4.0f*den_mid[Nx-3]+1.0f*den_mid[Nx-4])/ (dx*dx);
	}
}
*/

/*
void GrapheneFluid1D::BohmSource(string grid) {
	float dx3  = dx*dx*dx;
	if(grid=="mid"){
		for ( int i = 2; i <= Nx-3 ; i++ )
		{
			d3_den[i] = ( -0.5f*Den[i-2]+ Den[i-1]  -Den[i+1] + 0.5f*Den[i+2] ) / dx3;
			//d3_den_mid[i] = 2.0f*( -1.0f*Den[i-2]+ 3.0f*Den[i-1]  -3.0f*Den[i+1] + Den[i+2] ) / dx3;
		}
		//d3_den_mid[0]= 2.0f*(-5.0f*Den[0]	 +7.0f*Den[1]    +Den[2]	-3.0f*Den[3] )/ dx3;
		//d3_den_mid[1]= 2.0f*(-5.0f*Den[1]	 +7.0f*Den[2]    +Den[3]	-3.0f*Den[4] )/ dx3;
		//d3_den_mid[Nx-1]=2.0f*(5.0f*Den[Nx-1]-7.0f*Den[Nx-2] -Den[Nx-3]	+3.0f*Den[Nx-4] )/ dx3;
		//d3_den_mid[Nx-2]=2.0f*(5.0f*Den[Nx-2]-7.0f*Den[Nx-3] -Den[Nx-4]	+3.0f*Den[Nx-5] )/ dx3;
	}
	if(grid=="main"){
		for ( int i = 2; i <= Nx-3 ; i++ )
		{
			d3_den_mid[i] = 0.5f*( -1.0f*den_mid[i-2]+ 2.0f*den_mid[i-1]  -2.0f*den_mid[i+1] + den_mid[i+2] ) / dx3;
			//d3_den[i] =  2.0f*( -1.0f*den_mid[i-2]+ 3.0f*den_mid[i-1]  -3.0f*den_mid[i+1] + den_mid[i+2] ) / dx3;
		}
		d3_den[0]= (-2.5f*den_mid[0]	+9.0f*den_mid[1]-12.0f*den_mid[2]	+7.0f*den_mid[3]-1.5f*den_mid[4] )/ dx3;
		d3_den[0]= (-2.5f*den_mid[0]	+9.0f*den_mid[1]-12.0f*den_mid[2]	+7.0f*den_mid[3]-1.5f*den_mid[4] )/ dx3;
		d3_den[Nx-1]=(2.5f*den_mid[Nx-1]	-9.0f*den_mid[Nx-2]	+12.0f*den_mid[Nx-3]	-7.0f*den_mid[Nx-4]	+1.5f*den_mid[Nx-5] )/ dx3;
		d3_den[Nx-1]=(2.5f*den_mid[Nx-1]	-9.0f*den_mid[Nx-2]	+12.0f*den_mid[Nx-3]	-7.0f*den_mid[Nx-4]	+1.5f*den_mid[Nx-5] )/ dx3;
		//d3_den[0]= 2.0f*(-5.0f*den_mid[0]	 +7.0f*den_mid[1]    +den_mid[2]	-3.0f*den_mid[3] )/ dx3;
		//d3_den[1]= 2.0f*(-5.0f*den_mid[1]	 +7.0f*den_mid[2]    +den_mid[3]	-3.0f*den_mid[4] )/ dx3;
		//d3_den[Nx-1]=2.0f*(5.0f*den_mid[Nx-1]-7.0f*den_mid[Nx-2] -den_mid[Nx-3]	+3.0f*den_mid[Nx-4] )/ dx3;
		//d3_den[Nx-2]=2.0f*(5.0f*den_mid[Nx-2]-7.0f*den_mid[Nx-3] -den_mid[Nx-4]	+3.0f*den_mid[Nx-5] )/ dx3;
	}
}
*/