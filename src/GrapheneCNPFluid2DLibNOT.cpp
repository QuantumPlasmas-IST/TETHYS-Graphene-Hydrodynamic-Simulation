/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include <fftw3.h>
#include "includes/GrapheneCNPFluid2DLibNOT.h"
#include "external/polylogarithm-6.11.0/src/cpp/Li2.hpp"

using namespace polylogarithm;

GrapheneCNP2D::GrapheneCNP2D(SetUpParametersCNP &input_parameters) : FluidCNP2D(input_parameters) {
	//BLG
	tmp = input_parameters.InitialTemperature; //initial temperature
	vel_0 = input_parameters.InitialVelocity; //initial fluid velocity
	E0 = input_parameters.ExternalEfield; //external electric field
	charge_imb = input_parameters.ChargeImbalance; //charge imbalance

	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency;
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency;
	therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity
	char buffer [100];
	sprintf (buffer, "T=%.2fv0=%.2fE0=%.3fx=%.3fl=%.3fwc=%.2fvars=%d", tmp, vel_0, E0, charge_imb, col_freq, cyc_freq, num_vars);
	file_infix = buffer;
	vmax_x = 10.;
	vmax_y = 10.;
	time_counter = 0.;
}

void GrapheneCNP2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

void GrapheneCNP2D::MassFluxToVelocity(){
float den;
	for(int c=0; c <= Nx * Ny - 1; c++){
		den = Den[c];
		VelX[c] = FlxX[c] / sqrt(den*den*den);
		VelY[c] = FlxY[c] / sqrt(den*den*den);
	}
}

void GrapheneCNP2D::CflCondition(){ // Eventual redefinition
	dx = lengX / ( float ) (Nx - 1);
	dy = lengY / ( float ) (Ny - 1);
	dt = 0.5*min(dx/vmax_x, dy/vmax_y);
}



float GrapheneCNP2D::DensityFluxX(GridPoint p, char side) {
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

float GrapheneCNP2D::DensityFluxY(GridPoint p, char side) {
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

float GrapheneCNP2D::XMomentumFluxX(GridPoint p, char side) {
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

float GrapheneCNP2D::XMomentumFluxY(GridPoint p, char side) {
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


float GrapheneCNP2D::YMomentumFluxY(GridPoint p, char side) {
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

float GrapheneCNP2D::YMomentumFluxX(GridPoint p, char side) {
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


float GrapheneCNP2D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return 0.0f;
}

float GrapheneCNP2D::XMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_x  - cyc_freq*flx_y/sqrt(n);
}

float GrapheneCNP2D::YMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_y  + cyc_freq*flx_x/sqrt(n);
}


GrapheneCNP2D::~GrapheneCNP2D(){
delete[] Den;
delete[] Den_h;

delete[] VelX;
delete[] VelX_h;

delete[] VelY;
delete[] VelY_h;

delete[] FlxX;
delete[] FlxY;

delete[] CurX;
delete[] CurX_h;

delete[] CurY;
delete[] CurY_h;

delete[] E;
delete[] E_h;

delete[] P;
delete[] P_h;

delete[] Tmp;
delete[] Tmp_h;

delete[] cur_totX;
delete[] cur_totY;

delete[] EX;
delete[] EY;

delete[] den_E_energy;
delete[] den_E_th;
delete[] den_E_kin;

delete[] UL_x;
delete[] UR_x;

delete[] UL_y;
delete[] UR_y;

delete[] F;
delete[] G;

delete[] den_half;
delete[] den_h_half;

delete[] curX_half;
delete[] curY_half;

delete[] curX_h_half;
delete[] curY_h_half;

delete[] e_half;
delete[] e_h_half; 

delete[] p_half;
delete[] p_h_half;

delete[] t_half;
delete[] t_h_half;

delete[] phi;

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

float GrapheneCNP2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}

float
GrapheneCNP2D::TemperatureSource(float n, float flx_x, float flx_y, float den_grad_x, float den_grad_y, float mass, float s) {
	return vel_snd * vel_snd * (den_grad_x * flx_x / sqrt(n) + den_grad_y * flx_y / sqrt(n)  ) / (vel_fer * vel_fer);
}

//
// BLG METHODS
//

void GrapheneCNP2D::InitialCondRand2species(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();

	//Adimensional parameter equal to the ratio between chemical potential and temperature. Defines the charge imbalance
	float x = charge_imb;

	float r = log(1. + exp(x))/log(1. + exp(-x));
	float r_e = log(1. + exp(-x))/(log(1. + exp(x)) + log(1. + exp(-x)));

	//Initial values of temperature
	save_T_e = tmp;
	save_T_h = tmp;

	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;

			float noise_e = (float) rd()/maxrand;
			float noise_h = (float) rd()/maxrand;
			
			Tmp[m] = save_T_e;
        	Tmp_h[m] = save_T_h;

			Den[m] = Tmp[m]*log(1. + exp(x))/log(2.) + 0.005f*(noise_e - 0.5f);
			Den_h[m] = Tmp_h[m]*log(1. + exp(-x))/log(2.) + 0.005f*(noise_h - 0.5f);

			VelX[m] = -vel_0;
			VelY[m] = 0.f;
			
			VelX_h[m] = vel_0;
			VelY_h[m] = 0.f;

			CurX[m] = VelX[m]*Den[m];
			CurY[m] = VelY[m]*Den[m];
			CurX_h[m] = VelX_h[m]*Den_h[m];
			CurY_h[m] = VelY_h[m]*Den_h[m];

        	E[m] =  0.5*Den[m]*(VelX[m]*VelX[m] + VelY[m]*VelY[m]) - (Tmp[m]*Tmp[m]/log(2.))*Li2(1. - exp(log(2.)*Den[m]/Tmp[m]));
        	E_h[m] = 0.5*Den_h[m]*(VelX_h[m]*VelX_h[m] + VelY_h[m]*VelY_h[m]) - (Tmp_h[m]*Tmp_h[m]/log(2.))*Li2(1. - exp(log(2.)*Den_h[m]/Tmp_h[m]));

        	P[m] = -(Tmp[m]*Tmp[m]/log(2.))*Li2(1. - exp(log(2.)*Den[m]/Tmp[m]));
        	P_h[m] = -(Tmp_h[m]*Tmp_h[m]/log(2.))*Li2(1. - exp(log(2.)*Den_h[m]/Tmp_h[m]));
		}
		Den[i] = Den[i + (Ny - 2)*Nx];
		Den[i + (Ny - 1)*Nx] = Den[i + Nx]; 
		
		Den_h[i] = Den_h[i + (Ny - 2)*Nx];
		Den_h[i + (Ny - 1)*Nx] = Den_h[i + Nx];

		CurX[i] = CurX[i + (Ny - 2)*Nx];
		CurX[i + (Ny - 1)*Nx] = CurX[i + Nx]; 
		
		CurY[i] = CurY[i + (Ny - 2)*Nx];
		CurY[i + (Ny - 1)*Nx] = CurY[i + Nx]; 

		CurX_h[i] = CurX_h[i + (Ny - 2)*Nx];
		CurX_h[i + (Ny - 1)*Nx] = CurX_h[i + Nx]; 
		
		CurY_h[i] = CurY_h[i + (Ny - 2)*Nx];
		CurY_h[i + (Ny - 1)*Nx] = CurY_h[i + Nx];

		E[i] = E[i + (Ny - 2)*Nx];
		E[i + (Ny - 1)*Nx] = E[i + Nx]; 
		
		E_h[i] = E_h[i + (Ny - 2)*Nx];
		E_h[i + (Ny - 1)*Nx] = E_h[i + Nx];

		Tmp[i] = Tmp[i + (Ny - 2)*Nx];
		Tmp[i + (Ny - 1)*Nx] = Tmp[i + Nx]; 
		
		Tmp_h[i] = Tmp_h[i + (Ny - 2)*Nx];
		Tmp_h[i + (Ny - 1)*Nx] = Tmp_h[i + Nx];

		P[i] = P[i + (Ny - 2)*Nx];
		P[i + (Ny - 1)*Nx] = P[i + Nx]; 
		
		P_h[i] = P_h[i + (Ny - 2)*Nx];
		P_h[i + (Ny - 1)*Nx] = P_h[i + Nx];
	}
	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		
		Den[m0] = Den[m0 + Nx - 2];
		Den[mN] = Den[m0 + 1];
		
		Den_h[m0] = Den_h[m0 + Nx - 2];
		Den_h[mN] = Den_h[m0 + 1];

		CurX[m0] = CurX[m0 + Nx - 2];
		CurX[mN] = CurX[m0 + 1]; 
		
		CurY[m0] = CurY[m0 + Nx - 2];
		CurY[mN] = CurY[m0 + 1]; 

		CurX_h[m0] = CurX_h[m0 + Nx - 2];
		CurX_h[mN] = CurX_h[m0 + 1]; 
		
		CurY_h[m0] = CurY_h[m0 + Nx - 2];
		CurY_h[mN] = CurY_h[m0 + 1];

		E[m0] = E[m0 + Nx - 2];
		E[mN] = E[m0 + 1];
		
		E_h[m0] = E_h[m0 + Nx - 2];
		E_h[mN] = E_h[m0 + 1];

		Tmp[m0] = Tmp[m0 + Nx - 2];
		Tmp[mN] = Tmp[m0 + 1];
		
		Tmp_h[m0] = Tmp_h[m0 + Nx - 2];
		Tmp_h[mN] = Tmp_h[m0 + 1];

		P[m0] = P[m0 + Nx - 2];
		P[mN] = P[m0 + 1];
		
		P_h[m0] = P_h[m0 + Nx - 2];
		P_h[mN] = P_h[m0 + 1];
	}
}

float GrapheneCNP2D::Get_slope_x(const float* arr, int m, int flag){
	float forward_dif_x = arr[m + 1] - arr[m];
	float central_dif_x = 0.5*(arr[m + 1] - arr[m - 1]);
	float backward_dif_x = arr[m] - arr[m - 1];
	
	float r = forward_dif_x/backward_dif_x;
	float k = (float)flag*1.f/3.f;
	float beta = 1.;

	float delta_minus = backward_dif_x*max(0.f, min(1.f, beta*r));
	float delta_plus = backward_dif_x*max(0.f, min(beta, r));

	float slope = ((1.f - k)*delta_minus + (1.f + k)*delta_plus);


	return slope;
}

float GrapheneCNP2D::Get_slope_y(const float* arr, int m, int flag){
	float forward_dif_y = arr[m + Nx] - arr[m];
	float central_dif_y = 0.5*(arr[m + Nx] - arr[m - Nx]);
	float backward_dif_y = arr[m] - arr[m - Nx];
	
	float r = forward_dif_y/backward_dif_y;
	float k = (float)flag*1.f/3.f;
	float beta = 1.;

	float delta_minus = backward_dif_y*max(0.f, min(1.f, beta*r));
	float delta_plus = backward_dif_y*max(0.f, min(beta, r));

	float slope = ((1.f - k)*delta_minus + (1.f + k)*delta_plus);

	return slope;
}

void GrapheneCNP2D::Reconstruction(const float* n, const float* n_h, const float* cur_x, const float* cur_y, const float* cur_x_h, const float* cur_y_h,
                                   const float* e, const float* e_h, const float* p, const float* p_h, const float* t, const float* t_h){
	//Compute the slopes along x direction at each cell centroid using a midmod slope limiter and compute value of quantities
	//at the interface along x
	#pragma omp parallel for  default(none) shared(UL_x, UR_x, UL_y, UR_y, n, n_h, cur_x, cur_y, cur_x_h, cur_y_h, e, e_h, p, p_h, t, t_h)
	for (int j = Ny - 2; j >= 1; --j){
		for (int i = 1; i < Nx - 1; ++i){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			float slope_x_L, slope_x_R, slope_y_L, slope_y_R;

			//0
			slope_x_L = Get_slope_x(n, m, 1);
			slope_x_R = Get_slope_x(n, m, -1);
			slope_y_L = Get_slope_y(n, m, 1);
			slope_y_R = Get_slope_y(n, m, -1);
			
			UL_x[m_mid_x] = n[m] + 0.25*slope_x_L;
			UR_x[m_mid_x - 1] = n[m] - 0.25*slope_x_R;

			UL_y[m_mid_y] = n[m] + 0.25*slope_y_L;
			UR_y[m_mid_y - 1] = n[m] - 0.25*slope_y_R;

			//1
			slope_x_L = Get_slope_x(n_h, m, 1);
			slope_x_R = Get_slope_x(n_h, m, -1);
			slope_y_L = Get_slope_y(n_h, m, 1);
			slope_y_R = Get_slope_y(n_h, m, -1);

			UL_x[(Nx - 1)*(Ny - 2) + m_mid_x] = n_h[m] + 0.25*slope_x_L;
			UR_x[(Nx - 1)*(Ny - 2) + m_mid_x - 1] = n_h[m] - 0.25*slope_x_R;

			UL_y[(Ny - 1)*(Nx - 2) + m_mid_y] = n_h[m] + 0.25*slope_y_L;
			UR_y[(Ny - 1)*(Nx - 2) + m_mid_y - 1] = n_h[m] - 0.25*slope_y_R;

			//2
			slope_x_L = Get_slope_x(cur_x, m, 1);
			slope_x_R = Get_slope_x(cur_x, m, -1);
			slope_y_L = Get_slope_y(cur_x, m, 1);
			slope_y_R = Get_slope_y(cur_x, m, -1);

			UL_x[2*(Nx - 1)*(Ny - 2) + m_mid_x] = cur_x[m] + 0.25*slope_x_L;
			UR_x[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = cur_x[m] - 0.25*slope_x_R;

			UL_y[2*(Ny - 1)*(Nx - 2) + m_mid_y] = cur_x[m] + 0.25*slope_y_L;
			UR_y[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = cur_x[m] - 0.25*slope_y_R;

			//3
			slope_x_L = Get_slope_x(cur_y, m, 1);
			slope_x_R = Get_slope_x(cur_y, m, -1);
			slope_y_L = Get_slope_y(cur_y, m, 1);
			slope_y_R = Get_slope_y(cur_y, m, -1);

			UL_x[3*(Nx - 1)*(Ny - 2) + m_mid_x] = cur_y[m] + 0.25*slope_x_L;
			UR_x[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = cur_y[m] - 0.25*slope_x_R;

			UL_y[3*(Ny - 1)*(Nx - 2) + m_mid_y] = cur_y[m] + 0.25*slope_y_L;
			UR_y[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = cur_y[m] - 0.25*slope_y_R;

			//4
			slope_x_L = Get_slope_x(cur_x_h, m, 1);
			slope_x_R = Get_slope_x(cur_x_h, m, -1);
			slope_y_L = Get_slope_y(cur_x_h, m, 1);
			slope_y_R = Get_slope_y(cur_x_h, m, -1);

			UL_x[4*(Nx - 1)*(Ny - 2) + m_mid_x] = cur_x_h[m] + 0.25*slope_x_L;
			UR_x[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = cur_x_h[m] - 0.25*slope_x_R;

			UL_y[4*(Ny - 1)*(Nx - 2) + m_mid_y] = cur_x_h[m] + 0.25*slope_y_L;
			UR_y[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = cur_x_h[m] - 0.25*slope_y_R;

			//5
			slope_x_L = Get_slope_x(cur_y_h, m, 1);
			slope_x_R = Get_slope_x(cur_y_h, m, -1);
			slope_y_L = Get_slope_y(cur_y_h, m, 1);
			slope_y_R = Get_slope_y(cur_y_h, m, -1);

			UL_x[5*(Nx - 1)*(Ny - 2) + m_mid_x] = cur_y_h[m] + 0.25*slope_x_L;
			UR_x[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = cur_y_h[m] - 0.25*slope_x_R;

			UL_y[5*(Ny - 1)*(Nx - 2) + m_mid_y] = cur_y_h[m] + 0.25*slope_y_L;
			UR_y[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = cur_y_h[m] - 0.25*slope_y_R;

			if (num_vars > 6){
				//Electron pressure
				slope_x_L = Get_slope_x(p, m, 1);
				slope_x_R = Get_slope_x(p, m, -1);
				slope_y_L = Get_slope_y(p, m, 1);
				slope_y_R = Get_slope_y(p, m, -1);

				UL_x[8*(Nx - 1)*(Ny - 2) + m_mid_x] = p[m] + 0.25*slope_x_L;
				UR_x[8*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = p[m] - 0.25*slope_x_R;

				UL_y[8*(Ny - 1)*(Nx - 2) + m_mid_y] = p[m] + 0.25*slope_y_L;
				UR_y[8*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = p[m] - 0.25*slope_y_R;

				//Hole pressure
				slope_x_L = Get_slope_x(p_h, m, 1);
				slope_x_R = Get_slope_x(p_h, m, -1);
				slope_y_L = Get_slope_y(p_h, m, 1);
				slope_y_R = Get_slope_y(p_h, m, -1);

				UL_x[9*(Nx - 1)*(Ny - 2) + m_mid_x] = p_h[m] + 0.25*slope_x_L;
				UR_x[9*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = p_h[m] - 0.25*slope_x_R;

				UL_y[9*(Ny - 1)*(Nx - 2) + m_mid_y] = p_h[m] + 0.25*slope_y_L;
				UR_y[9*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = p_h[m] - 0.25*slope_y_R;

				//Compute Electron energy densities from pressure
				UL_x[6*(Nx - 1)*(Ny - 2) + m_mid_x] = UL_x[8*(Nx - 1)*(Ny - 2) + m_mid_x] 
				+ 0.5*(UL_x[2*(Nx - 1)*(Ny - 2) + m_mid_x]*UL_x[2*(Nx - 1)*(Ny - 2) + m_mid_x] + 
				UL_x[3*(Nx - 1)*(Ny - 2) + m_mid_x]*UL_x[3*(Nx - 1)*(Ny - 2) + m_mid_x])/UL_x[m_mid_x];
				
				UR_x[6*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = UR_x[8*(Nx - 1)*(Ny - 2) + m_mid_x - 1]
				 + 0.5*(UR_x[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1]*UR_x[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1] + 
				UR_x[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1]*UR_x[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1])/UR_x[m_mid_x - 1];

				UL_y[6*(Ny - 1)*(Nx - 2) + m_mid_y] = UL_y[8*(Ny - 1)*(Nx - 2) + m_mid_y] 
				+ 0.5*(UL_y[2*(Ny - 1)*(Nx - 2) + m_mid_y]*UL_y[2*(Ny - 1)*(Nx - 2) + m_mid_y] + 
				UL_y[3*(Ny - 1)*(Nx - 2) + m_mid_y]*UL_y[3*(Ny - 1)*(Nx - 2) + m_mid_y])/UL_y[m_mid_y];

				UR_y[6*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = UR_y[8*(Ny - 1)*(Nx - 2) + m_mid_y - 1] 
				+ 0.5*(UR_y[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1]*UR_y[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1] + 
				UR_y[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1]*UR_y[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1])/UR_y[m_mid_y - 1];

				//Compute Hole energy densities from pressure
				UL_x[7*(Nx - 1)*(Ny - 2) + m_mid_x] = UL_x[9*(Nx - 1)*(Ny - 2) + m_mid_x] 
				+ 0.5*(UL_x[4*(Nx - 1)*(Ny - 2) + m_mid_x]*UL_x[4*(Nx - 1)*(Ny - 2) + m_mid_x] + 
				UL_x[5*(Nx - 1)*(Ny - 2) + m_mid_x]*UL_x[5*(Nx - 1)*(Ny - 2) + m_mid_x])/UL_x[(Nx - 1)*(Ny - 2) + m_mid_x];
				
				UR_x[7*(Nx - 1)*(Ny - 2) + m_mid_x - 1] = UR_x[9*(Nx - 1)*(Ny - 2) + m_mid_x - 1]
				 + 0.5*(UR_x[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1]*UR_x[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1] + 
				UR_x[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1]*UR_x[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1])/UR_x[(Nx - 1)*(Ny - 2) + m_mid_x - 1];

				UL_y[7*(Ny - 1)*(Nx - 2) + m_mid_y] = UL_y[9*(Ny - 1)*(Nx - 2) + m_mid_y] 
				+ 0.5*(UL_y[4*(Ny - 1)*(Nx - 2) + m_mid_y]*UL_y[4*(Ny - 1)*(Nx - 2) + m_mid_y] + 
				UL_y[5*(Ny - 1)*(Nx - 2) + m_mid_y]*UL_y[5*(Ny - 1)*(Nx - 2) + m_mid_y])/UL_y[(Ny - 1)*(Nx - 2) + m_mid_y];

				UR_y[7*(Ny - 1)*(Nx - 2) + m_mid_y - 1] = UR_y[9*(Ny - 1)*(Nx - 2) + m_mid_y - 1] 
				+ 0.5*(UR_y[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1]*UR_y[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1] + 
				UR_y[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1]*UR_y[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1])/UR_y[(Ny - 1)*(Nx - 2) + m_mid_y - 1];
			}
			
			for (int flag = 0; flag < num_vars; ++flag){
				UL_y[flag*(Ny - 1)*(Nx - 2) + (i - 1)*(Ny - 1)] = UL_y[flag*(Ny - 1)*(Nx - 2) + (i - 1)*(Ny - 1) + Ny - 2];
				UR_y[flag*(Ny - 1)*(Nx - 2) + (i - 1)*(Ny - 1) + Ny - 2] = UR_y[flag*(Ny - 1)*(Nx - 2) + (i - 1)*(Ny - 1)];
			}
		}
		for (int flag = 0; flag < num_vars; ++flag){
			UL_x[flag*(Nx - 1)*(Ny - 2) + (j - 1)*(Nx - 1)] = UL_x[flag*(Nx - 1)*(Ny - 2) + (j - 1)*(Nx - 1) + Nx - 2];
			UR_x[flag*(Nx - 1)*(Ny - 2) + (j - 1)*(Nx - 1) + Nx - 2] = UR_x[flag*(Nx - 1)*(Ny - 2) + (j - 1)*(Nx - 1)];
		}
	}
}

float GrapheneCNP2D::HLLFlux(float FL, float FR, float UL, float UR, float UL_star, float UR_star, float SL, float SR, float SM){
	float flux;

	//Intermediate fluxes
	float FL_star = SL*(UL_star - UL) + FL;
	float FR_star = SR*(UR_star - UR) + FR;

	float F_star = 0.5*(FL + FR) + 0.5*(SL*(UL_star - UL) + fabs(SM)*(UL_star - UR_star) + SR*(UR_star - UR));

	if (SL > 0.f){
		flux = FL;
	}
	else if (SL <= 0.f && SM > 0.f){
		flux = FL_star;
	}
	else if (SM <= 0.f && SR > 0.f){
		flux = FR_star;
	}
	else if (SR <= 0.f){
		flux = FR;
	}

	return flux;
}

float GrapheneCNP2D::GetTmp(float p, float den, float trial){
	float aux = 144.*p*p - 72.*den*den*p*log(2.) - 7.*den*den*den*den*log(2)*log(2.);
	float x = -(3.*(-12.*p + sqrt(144.*p*p - 72.*den*den*p*log(2.) - 7.*den*den*den*den*log(2)*log(2.)) + den*den*log(8.)))/(2.*den*den*log(2.)*log(2.));
	float tmp;

	//If the linear approximation yields a physical result, we use it. If not, we solve the non-linear equation numerically
	if (aux >= 0. && x > 0. && p > 0.){
		tmp = den/x; 
	}
	else{
		float tol = 1e-5;
		float tmp = tol;
			
		//Initial trial solution to initialize the method
		float t0 = trial;

		//Iterate until the specified tolerance is achieved or after a maximum number of iterations is reached
		int max_count = 10000;
		float err;

		if (p == 0.){
			cout << "p" << "\t" << p << endl;
		}

		for (int i = 0; i < max_count; ++i){
			float t = sqrt(-p*log(2.)/Li2(1. - exp(log(2.)*den/t0)));
			err = fabs(t - t0);

			if (err < tol && t > 0.){
				tmp = t;
				break;
			}

			if (t < tol){
				tmp = t;
				break;
			}
			t0 = t;
		}
		if (isnan(t0)){
			if (p <= 0.)
				cout << "p" << "\t" << p << endl;
			tmp = tol;
		}
	}

	return tmp;
}

void GrapheneCNP2D::ComputeFluxX(){
	//Compute the flux at the interfaces along x direction
	#pragma omp parallel for default(none) shared(F, UL_x, UR_x)
	for (int m = 0; m < (Ny - 2)*(Nx - 1); ++m){
		//Electron left and right primitive variable states
		float cs_e_L = sqrt(log(2.)*UL_x[m]/(1. - exp(-log(2.)*UL_x[m]/save_T_e)));
		float cs_e_R = sqrt(log(2.)*UR_x[m]/(1. - exp(-log(2.)*UR_x[m]/save_T_e)));
		
		float vx_e_L = UL_x[2*(Nx - 1)*(Ny - 2) + m]/UL_x[m];
		float vx_e_R = UR_x[2*(Nx - 1)*(Ny - 2) + m]/UR_x[m];

		float vy_e_L = UL_x[3*(Nx - 1)*(Ny - 2) + m]/UL_x[m];
		float vy_e_R = UR_x[3*(Nx - 1)*(Ny - 2) + m]/UR_x[m];

		float p_e_L = -(save_T_e*save_T_e/log(2.))*Li2(1. - exp(log(2.)*UL_x[m]/save_T_e));
		float p_e_R = -(save_T_e*save_T_e/log(2.))*Li2(1. - exp(log(2.)*UR_x[m]/save_T_e));		

		//Hole left and right primitive variable states
		float cs_h_L = sqrt(log(2.)*UL_x[(Nx - 1)*(Ny - 2) + m]/(1. - exp(-log(2.)*UL_x[(Nx - 1)*(Ny - 2) + m]/save_T_h)));
		float cs_h_R = sqrt(log(2.)*UR_x[(Nx - 1)*(Ny - 2) + m]/(1. - exp(-log(2.)*UR_x[(Nx - 1)*(Ny - 2) + m]/save_T_h)));
		
		float vx_h_L = UL_x[4*(Nx - 1)*(Ny - 2) + m]/UL_x[(Nx - 1)*(Ny - 2) + m];
		float vx_h_R = UR_x[4*(Nx - 1)*(Ny - 2) + m]/UR_x[(Nx - 1)*(Ny - 2) + m];
		
		float vy_h_L = UL_x[5*(Nx - 1)*(Ny - 2) + m]/UL_x[(Nx - 1)*(Ny - 2) + m];
		float vy_h_R = UR_x[5*(Nx - 1)*(Ny - 2) + m]/UR_x[(Nx - 1)*(Ny - 2) + m];

		float p_h_L = -(save_T_h*save_T_h/log(2.))*Li2(1. - exp(log(2.)*UL_x[(Nx - 1)*(Ny - 2) + m]/save_T_h));
		float p_h_R = -(save_T_h*save_T_h/log(2.))*Li2(1. - exp(log(2.)*UR_x[(Nx - 1)*(Ny - 2) + m]/save_T_h));

		float t_e_L, t_e_R, t_h_L, t_h_R;

		if (num_vars > 6){
			p_e_L = UL_x[8*(Nx - 1)*(Ny - 2) + m];
			p_e_R = UR_x[8*(Nx - 1)*(Ny - 2) + m];
			p_h_L = UL_x[9*(Nx - 1)*(Ny - 2) + m];
			p_h_R = UR_x[9*(Nx - 1)*(Ny - 2) + m];

			float t1 = GetTmp(p_e_L, UL_x[m], 1.0f);
			float t2 = GetTmp(p_e_R, UR_x[m], 1.0f);
			float t3 = GetTmp(p_h_L, UL_x[(Nx - 1)*(Ny - 2) + m], 1.0f);
			float t4 = GetTmp(p_h_R, UR_x[(Nx - 1)*(Ny - 2) + m], 1.0f);

			t_e_L = t1;
			t_e_R = t2;
			t_h_L = t3;
			t_h_R = t4;

			cs_e_L = sqrt(log(2.)*UL_x[m]/(1. - exp(-log(2.)*UL_x[m]/t_e_L)));
			cs_e_R = sqrt(log(2.)*UR_x[m]/(1. - exp(-log(2.)*UR_x[m]/t_e_R)));
			cs_h_L = sqrt(log(2.)*UL_x[(Nx - 1)*(Ny - 2) + m]/(1. - exp(-log(2.)*UL_x[(Nx - 1)*(Ny - 2) + m]/t_h_L)));
			cs_h_R = sqrt(log(2.)*UR_x[(Nx - 1)*(Ny - 2) + m]/(1. - exp(-log(2.)*UR_x[(Nx - 1)*(Ny - 2) + m]/t_h_R)));	
		}

		//Roe average velocities
		float u_e = (vx_e_L*sqrt(UL_x[m]) + vx_e_R*sqrt(UR_x[m]))/(sqrt(UL_x[m]) + sqrt(UR_x[m]));
		
		float c_e = sqrt((cs_e_L*cs_e_L*sqrt(UL_x[m]) + cs_e_R*cs_e_R*sqrt(UR_x[m]))/(sqrt(UL_x[m]) + sqrt(UR_x[m])) + 
		0.5*sqrt(UL_x[m]*UR_x[m])*(vx_e_L - vx_e_R)*(vx_e_L - vx_e_R)/((sqrt(UL_x[m]) + sqrt(UR_x[m]))*(sqrt(UL_x[m]) + sqrt(UR_x[m]))));
		
		float u_h = (vx_h_L*sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) + vx_h_R*sqrt(UR_x[(Nx - 1)*(Ny - 2) + m]))/(sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) 
		+ sqrt(UR_x[(Nx - 1)*(Ny - 2) + m]));
		
		float c_h = sqrt((cs_h_L*cs_h_L*sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) + cs_h_R*cs_h_R*sqrt(UR_x[(Nx - 1)*(Ny - 2) + m]))/(sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) 
		+ sqrt(UR_x[(Nx - 1)*(Ny - 2) + m])) + 0.5*sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]*UR_x[(Nx - 1)*(Ny - 2) + m])
		*(vx_h_L - vx_h_R)*(vx_h_L - vx_h_R)/((sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) + sqrt(UR_x[(Nx - 1)*(Ny - 2) + m]))*(sqrt(UL_x[(Nx - 1)*(Ny - 2) + m]) 
		+ sqrt(UR_x[(Nx - 1)*(Ny - 2) + m]))));
		
		//Estimate the characteristic wave velocities
		float SL_e = min(vx_e_L - cs_e_L, u_e - c_e);
		float SR_e = max(vx_e_R + cs_e_R, u_e + c_e);
		float SL_h = min(vx_h_L - cs_h_L, u_h - c_h);
		float SR_h = max(vx_h_R + cs_h_R, u_h + c_h);

		//Characteristic wave velocity of the intermediate state
		float SM_e = ((SR_e - vx_e_R)*UR_x[m]*vx_e_R - (SL_e - vx_e_L)*UL_x[m]*vx_e_L - p_e_R + p_e_L)/((SR_e - vx_e_R)*UR_x[m] - (SL_e - vx_e_L)*UL_x[m]);
		float SM_h = ((SR_h - vx_h_R)*UR_x[(Nx - 1)*(Ny - 2) + m]*vx_h_R - (SL_h - vx_h_L)*UL_x[(Nx - 1)*(Ny - 2) + m]*vx_h_L 
		- p_h_R + p_h_L)/((SR_h - vx_h_R)*UR_x[(Nx - 1)*(Ny - 2) + m] - (SL_h - vx_h_L)*UL_x[(Nx - 1)*(Ny - 2) + m]); 

		//Density of the intermediate state
		float n_starL_e = UL_x[m]*(SL_e - vx_e_L)/(SL_e - SM_e);
		float n_starR_e = UR_x[m]*(SR_e - vx_e_R)/(SR_e - SM_e);
		float n_starL_h = UL_x[(Nx - 1)*(Ny - 2) + m]*(SL_h - vx_h_L)/(SL_h - SM_h);
		float n_starR_h = UR_x[(Nx - 1)*(Ny - 2) + m]*(SR_h - vx_h_R)/(SR_h - SM_h);

		//CurX of the intermediate state
		float CurX_starL_e = n_starL_e*SM_e;
		float CurX_starR_e = n_starR_e*SM_e;
		float CurX_starL_h = n_starL_h*SM_h;
		float CurX_starR_h = n_starR_h*SM_h;

		//CurY of the intermediate state
		float CurY_starL_e = n_starL_e*vy_e_L;
		float CurY_starR_e = n_starR_e*vy_e_R;
		float CurY_starL_h = n_starL_h*vy_h_L;
		float CurY_starR_h = n_starR_h*vy_h_R;

		float eig_max_L = max(fabs(vx_e_L) + cs_e_L, fabs(vx_h_L) + cs_h_L);
		float eig_max_R = max(fabs(vx_e_R) + cs_e_R, fabs(vx_h_R) + cs_h_R);
		float eig_max = max(eig_max_L, eig_max_R);

		if (vmax_x < eig_max){
			vmax_x = eig_max;
		}

		//Electron density
		F[m] = HLLFlux(DensityFluxX2species(UL_x[2*(Nx - 1)*(Ny - 2) + m]), DensityFluxX2species(UR_x[2*(Nx - 1)*(Ny - 2) + m]), 
			UL_x[m], UR_x[m], n_starL_e, n_starR_e, SL_e, SR_e, SM_e);
		
		//Hole density
		F[(Nx - 1)*(Ny - 2) + m] = HLLFlux(DensityFluxX2species(UL_x[4*(Nx - 1)*(Ny - 2) + m]), DensityFluxX2species(UR_x[4*(Nx - 1)*(Ny - 2) + m]), 
			UL_x[(Nx - 1)*(Ny - 2) + m], UR_x[(Nx - 1)*(Ny - 2) + m], n_starL_h, n_starR_h, SL_h, SR_h, SM_h);
		
		//Electron current x component of the Xflux
		F[2*(Nx - 1)*(Ny - 2) + m] = HLLFlux(XCurrentFluxX2species(UL_x[m], UL_x[2*(Nx - 1)*(Ny - 2) + m], p_e_L), 
			XCurrentFluxX2species(UR_x[m], UR_x[2*(Nx - 1)*(Ny - 2) + m], p_e_R), UL_x[2*(Nx - 1)*(Ny - 2) + m], UR_x[2*(Nx - 1)*(Ny - 2) + m], 
			CurX_starL_e, CurX_starR_e, SL_e, SR_e, SM_e);
		
		//Electron current y component of the Xflux
		F[3*(Nx - 1)*(Ny - 2) + m] = HLLFlux(XCurrentFluxY2species(UL_x[m], UL_x[2*(Nx - 1)*(Ny - 2) + m], UL_x[3*(Nx - 1)*(Ny - 2) + m]),
			XCurrentFluxY2species(UR_x[m], UR_x[2*(Nx - 1)*(Ny - 2) + m], UR_x[3*(Nx - 1)*(Ny - 2) + m]), 
			UL_x[3*(Nx - 1)*(Ny - 2) + m], UR_x[3*(Nx - 1)*(Ny - 2) + m], CurY_starL_e, CurY_starR_e, SL_e, SR_e, SM_e);

		//Hole current x component of the Xflux
		F[4*(Nx - 1)*(Ny - 2) + m] = HLLFlux(XCurrentFluxX2species(UL_x[(Nx - 1)*(Ny - 2) + m], UL_x[4*(Nx - 1)*(Ny - 2) + m], p_h_L), 
			XCurrentFluxX2species(UR_x[(Nx - 1)*(Ny - 2) + m], UR_x[4*(Nx - 1)*(Ny - 2) + m], p_h_R), UL_x[4*(Nx - 1)*(Ny - 2) + m], UR_x[4*(Nx - 1)*(Ny - 2) + m], 
			CurX_starL_h, CurX_starR_h, SL_h, SR_h, SM_h);
		
		//Hole current y component of the Xflux
		F[5*(Nx - 1)*(Ny - 2) + m] = HLLFlux(XCurrentFluxY2species(UL_x[(Nx - 1)*(Ny - 2) + m], UL_x[4*(Nx - 1)*(Ny - 2) + m], UL_x[5*(Nx - 1)*(Ny - 2) + m]), 
			XCurrentFluxY2species(UR_x[(Nx - 1)*(Ny - 2) + m], UR_x[4*(Nx - 1)*(Ny - 2) + m], UR_x[5*(Nx - 1)*(Ny - 2) + m]), 
			UL_x[5*(Nx - 1)*(Ny - 2) + m], UR_x[5*(Nx - 1)*(Ny - 2) + m], CurY_starL_h, CurY_starR_h, SL_h, SR_h, SM_h);

		if (num_vars > 6){
			//Energy density of intermediate state
			float E_starL_e = ((SL_e - vx_e_L)*UL_x[6*(Nx - 1)*(Ny - 2) + m] + p_e_L*(SM_e - vx_e_L) + SM_e*UL_x[m]*(SM_e - vx_e_L)*(SL_e - vx_e_L))/(SL_e - SM_e);
			float E_starR_e = ((SR_e - vx_e_R)*UR_x[6*(Nx - 1)*(Ny - 2) + m] + p_e_R*(SM_e - vx_e_R) + SM_e*UR_x[m]*(SM_e - vx_e_R)*(SR_e - vx_e_R))/(SR_e - SM_e);
			float E_starL_h = ((SL_h - vx_h_L)*UL_x[7*(Nx - 1)*(Ny - 2) + m] + p_h_L*(SM_h - vx_h_L) + SM_h*UL_x[(Nx - 1)*(Ny - 2) + m]*(SM_h - vx_h_L)*(SL_h - vx_h_L))/(SL_h - SM_h);
			float E_starR_h = ((SR_h - vx_h_R)*UR_x[7*(Nx - 1)*(Ny - 2) + m] + p_h_R*(SM_h - vx_h_R) + SM_h*UR_x[(Nx - 1)*(Ny - 2) + m]*(SM_h - vx_h_R)*(SR_h - vx_h_R))/(SR_h - SM_h);

			//Electron internal energy
			F[6*(Nx - 1)*(Ny - 2) + m] = HLLFlux(EnergyFluxX2species(UL_x[m], UL_x[2*(Nx - 1)*(Ny - 2) + m], UL_x[6*(Nx - 1)*(Ny - 2) + m], p_e_L), 
				EnergyFluxX2species(UR_x[m], UR_x[2*(Nx - 1)*(Ny - 2) + m], UR_x[6*(Nx - 1)*(Ny - 2) + m], p_e_R), 
				UL_x[6*(Nx - 1)*(Ny - 2) + m], UR_x[6*(Nx - 1)*(Ny - 2) + m], E_starL_e, E_starR_e, SL_e, SR_e, SM_e);

			//Hole internal energy
			F[7*(Nx - 1)*(Ny - 2) + m] = HLLFlux(EnergyFluxX2species(UL_x[(Nx - 1)*(Ny - 2) + m], UL_x[4*(Nx - 1)*(Ny - 2) + m], UL_x[7*(Nx - 1)*(Ny - 2) + m], p_h_L), 
				EnergyFluxX2species(UR_x[(Nx - 1)*(Ny - 2) + m], UR_x[4*(Nx - 1)*(Ny - 2) + m], UR_x[7*(Nx - 1)*(Ny - 2) + m], p_h_R), 
				UL_x[7*(Nx - 1)*(Ny - 2) + m], UR_x[7*(Nx - 1)*(Ny - 2) + m], E_starL_h, E_starR_h, SL_h, SR_h, SM_h);
		}
	}
}

void GrapheneCNP2D::ComputeFluxY(){
	//Compute the flux at the interfaces along y direction
	#pragma omp parallel for  default(none) shared(G, UL_y, UR_y)
	for (int m = 0; m < (Nx - 2)*(Ny - 1); ++m){
		//Electron left and right primitive variable states
		float cs_e_L = sqrt(log(2.)*UL_y[m]/(1. - exp(-log(2.)*UL_y[m]/save_T_e)));
		float cs_e_R = sqrt(log(2.)*UR_y[m]/(1. - exp(-log(2.)*UR_y[m]/save_T_e)));
		
		float vy_e_L = UL_y[3*(Ny - 1)*(Nx - 2) + m]/UL_y[m];
		float vy_e_R = UR_y[3*(Ny - 1)*(Nx - 2) + m]/UR_y[m];

		float vx_e_L = UL_y[2*(Ny - 1)*(Nx - 2) + m]/UL_y[m];
		float vx_e_R = UR_y[2*(Ny - 1)*(Nx - 2) + m]/UR_y[m];

		float p_e_L = -(save_T_e*save_T_e/log(2.))*Li2(1. - exp(log(2.)*UL_y[m]/save_T_e));
		float p_e_R = -(save_T_e*save_T_e/log(2.))*Li2(1. - exp(log(2.)*UR_y[m]/save_T_e));

		//Hole left and right primitive variable states
		float cs_h_L = sqrt(log(2.)*UL_y[(Ny - 1)*(Nx - 2) + m]/(1. - exp(-log(2.)*UL_y[(Ny - 1)*(Nx - 2) + m]/save_T_h)));
		float cs_h_R = sqrt(log(2.)*UR_y[(Ny - 1)*(Nx - 2) + m]/(1. - exp(-log(2.)*UR_y[(Ny - 1)*(Nx - 2) + m]/save_T_h)));
		
		float vy_h_L = UL_y[5*(Ny - 1)*(Nx - 2) + m]/UL_y[(Ny - 1)*(Nx - 2) + m];
		float vy_h_R = UR_y[5*(Ny - 1)*(Nx - 2) + m]/UR_y[(Ny - 1)*(Nx - 2) + m];

		float vx_h_L = UL_y[4*(Ny - 1)*(Nx - 2) + m]/UL_y[(Ny - 1)*(Nx - 2) + m];
		float vx_h_R = UR_y[4*(Ny - 1)*(Nx - 2) + m]/UR_y[(Ny - 1)*(Nx - 2) + m];

		float p_h_L = -(save_T_h*save_T_h/log(2.))*Li2(1. - exp(log(2.)*UL_y[(Ny - 1)*(Nx - 2) + m]/save_T_h));
		float p_h_R = -(save_T_h*save_T_h/log(2.))*Li2(1. - exp(log(2.)*UR_y[(Ny - 1)*(Nx - 2) + m]/save_T_h));

		float t_e_L, t_e_R, t_h_L, t_h_R;

		if (num_vars > 6){
			p_e_L = UL_y[8*(Ny - 1)*(Nx - 2) + m];
			p_e_R = UR_y[8*(Ny - 1)*(Nx - 2) + m];
			p_h_L = UL_y[9*(Ny - 1)*(Nx - 2) + m];
			p_h_R = UR_y[9*(Ny - 1)*(Nx - 2) + m];

			float t1 = GetTmp(p_e_L, UL_y[m], 1.0f);
			float t2 = GetTmp(p_e_R, UR_y[m], 1.0f);
			float t3 = GetTmp(p_h_L, UL_y[(Ny - 1)*(Nx - 2) + m], 1.0f);
			float t4 = GetTmp(p_h_R, UR_y[(Ny - 1)*(Nx - 2) + m], 1.0f);

			t_e_L = t1;
			t_e_R = t2;
			t_h_L = t3;
			t_h_R = t4;

			cs_e_L = sqrt(log(2.)*UL_y[m]/(1. - exp(-log(2.)*UL_y[m]/t_e_L)));
			cs_e_R = sqrt(log(2.)*UR_y[m]/(1. - exp(-log(2.)*UR_y[m]/t_e_R)));
			cs_h_L = sqrt(log(2.)*UL_y[(Ny - 1)*(Nx - 2) + m]/(1. - exp(-log(2.)*UL_y[(Ny - 1)*(Nx - 2) + m]/t_h_L)));
			cs_h_R = sqrt(log(2.)*UR_y[(Ny - 1)*(Nx - 2) + m]/(1. - exp(-log(2.)*UR_y[(Ny - 1)*(Nx - 2) + m]/t_h_R)));	
		}

		//Roe average velocities
		float u_e = (vy_e_L*sqrt(UL_y[m]) + vy_e_R*sqrt(UR_y[m]))/(sqrt(UL_y[m]) + sqrt(UR_y[m]));
		
		float c_e = sqrt((cs_e_L*cs_e_L*sqrt(UL_y[m]) + cs_e_R*cs_e_R*sqrt(UR_y[m]))/(sqrt(UL_y[m]) + sqrt(UR_y[m])) + 
		0.5*sqrt(UL_y[m]*UR_y[m])*(vy_e_L - vy_e_R)*(vy_e_L - vy_e_R)/((sqrt(UL_y[m]) + sqrt(UR_y[m]))*(sqrt(UL_y[m]) + sqrt(UR_y[m]))));
		
		float u_h = (vy_h_L*sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) + vy_h_R*sqrt(UR_y[(Ny - 1)*(Nx - 2) + m]))/(sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) 
		+ sqrt(UR_y[(Ny - 1)*(Nx - 2) + m]));
		
		float c_h = sqrt((cs_h_L*cs_h_L*sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) + cs_h_R*cs_h_R*sqrt(UR_y[(Ny - 1)*(Nx - 2) + m]))/(sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) 
		+ sqrt(UR_y[(Ny - 1)*(Nx - 2) + m])) + 0.5*sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]*UR_y[(Ny - 1)*(Nx - 2) + m])
		*(vy_h_L - vy_h_R)*(vy_h_L - vy_h_R)/((sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) + sqrt(UR_y[(Ny - 1)*(Nx - 2) + m]))*(sqrt(UL_y[(Ny - 1)*(Nx - 2) + m]) 
		+ sqrt(UR_y[(Ny - 1)*(Nx - 2) + m]))));
		
		//Estimate the characteristic wave velocities
		float SL_e = min(vy_e_L - cs_e_L, u_e - c_e);
		float SR_e = max(vy_e_R + cs_e_R, u_e + c_e);
		float SL_h = min(vy_h_L - cs_h_L, u_h - c_h);
		float SR_h = max(vy_h_R + cs_h_R, u_h + c_h);

		//Characteristic wave velocity of the intermediate state
		float SM_e = ((SR_e - vy_e_R)*UR_y[m]*vy_e_R - (SL_e - vy_e_L)*UL_y[m]*vy_e_L - p_e_R + p_e_L)/((SR_e - vy_e_R)*UR_y[m] - (SL_e - vy_e_L)*UL_y[m]);
		float SM_h = ((SR_h - vy_h_R)*UR_y[(Ny - 1)*(Nx - 2) + m]*vy_h_R - (SL_h - vy_h_L)*UL_y[(Ny - 1)*(Nx - 2) + m]*vy_h_L 
		- p_h_R + p_h_L)/((SR_h - vy_h_R)*UR_y[(Ny - 1)*(Nx - 2) + m] - (SL_h - vy_h_L)*UL_y[(Ny - 1)*(Nx - 2) + m]); 

		//Density of the intermediate state
		float n_starL_e = UL_y[m]*(SL_e - vy_e_L)/(SL_e - SM_e);
		float n_starR_e = UR_y[m]*(SR_e - vy_e_R)/(SR_e - SM_e);
		float n_starL_h = UL_y[(Ny - 1)*(Nx - 2) + m]*(SL_h - vy_h_L)/(SL_h - SM_h);
		float n_starR_h = UR_y[(Ny - 1)*(Nx - 2) + m]*(SR_h - vy_h_R)/(SR_h - SM_h);

		//CurX of the intermediate state
		float CurX_starL_e = n_starL_e*vx_e_L;
		float CurX_starR_e = n_starR_e*vx_e_R;
		float CurX_starL_h = n_starL_h*vx_h_L;
		float CurX_starR_h = n_starR_h*vx_h_R;

		//CurY of the intermediate state
		float CurY_starL_e = n_starL_e*SM_e;
		float CurY_starR_e = n_starR_e*SM_e;
		float CurY_starL_h = n_starL_h*SM_h;
		float CurY_starR_h = n_starR_h*SM_h;

		float eig_max_L = max(fabs(vy_e_L) + cs_e_L, fabs(vy_h_L) + cs_h_L);
		float eig_max_R = max(fabs(vy_e_R) + cs_e_R, fabs(vy_h_R) + cs_h_R);
		float eig_max = max(eig_max_L, eig_max_R);

		if (vmax_y < eig_max){
			vmax_y = eig_max;
		}

		//Electron density
		G[m] = HLLFlux(DensityFluxY2species(UL_y[3*(Ny - 1)*(Nx - 2) + m]), DensityFluxY2species(UR_y[3*(Ny - 1)*(Nx - 2) + m]), 
			UL_y[m], UR_y[m], n_starL_e, n_starR_e, SL_e, SR_e, SM_e);
		
		//Hole density
		G[(Ny - 1)*(Nx - 2) + m] = HLLFlux(DensityFluxY2species(UL_y[5*(Ny - 1)*(Nx - 2) + m]), DensityFluxY2species(UR_y[5*(Ny - 1)*(Nx - 2) + m]),
			UL_y[(Ny - 1)*(Nx - 2) + m], UR_y[(Ny - 1)*(Nx - 2) + m], n_starL_h, n_starR_h, SL_h, SR_h, SM_h);
		
		//Electron x component of the Yflux
		G[2*(Ny - 1)*(Nx - 2) + m] = HLLFlux(YCurrentFluxX2species(UL_y[m], UL_y[2*(Ny - 1)*(Nx - 2) + m], UL_y[3*(Ny - 1)*(Nx - 2) + m]),
			YCurrentFluxX2species(UR_y[m], UR_y[2*(Ny - 1)*(Nx - 2) + m], UR_y[3*(Ny - 1)*(Nx - 2) + m]), 
			UL_y[2*(Ny - 1)*(Nx - 2) + m], UR_y[2*(Ny - 1)*(Nx - 2) + m], CurX_starL_e, CurX_starR_e, SL_e, SR_e, SM_e);
		
		//Electron y component of the Yflux
		G[3*(Ny - 1)*(Nx - 2) + m] = HLLFlux(YCurrentFluxY2species(UL_y[m], UL_y[3*(Ny - 1)*(Nx - 2) + m], p_e_L), 
			YCurrentFluxY2species(UR_y[m], UR_y[3*(Ny - 1)*(Nx - 2) + m], p_e_R), UL_y[3*(Ny - 1)*(Nx - 2) + m], UR_y[3*(Ny - 1)*(Nx - 2) + m], 
			CurY_starL_e, CurY_starR_e, SL_e, SR_e, SM_e);

		//Hole x component of the Yflux
		G[4*(Ny - 1)*(Nx - 2) + m] = HLLFlux(YCurrentFluxX2species(UL_y[(Ny - 1)*(Nx - 2) + m], UL_y[4*(Ny - 1)*(Nx - 2) + m], UL_y[5*(Ny - 1)*(Nx - 2) + m]), 
			YCurrentFluxX2species(UR_y[(Ny - 1)*(Nx - 2) + m], UR_y[4*(Ny - 1)*(Nx - 2) + m], UR_y[5*(Ny - 1)*(Nx - 2) + m]), 
			UL_y[4*(Ny - 1)*(Nx - 2) + m], UR_y[4*(Ny - 1)*(Nx - 2) + m], CurX_starL_h, CurX_starR_h, SL_h, SR_h, SM_h);
		
		//Hole y component of the Yflux
		G[5*(Ny - 1)*(Nx - 2) + m] = HLLFlux(YCurrentFluxY2species(UL_y[(Ny - 1)*(Nx - 2) + m], UL_y[5*(Ny - 1)*(Nx - 2) + m], p_h_L), 
			YCurrentFluxY2species(UR_y[(Ny - 1)*(Nx - 2) + m], UR_y[5*(Ny - 1)*(Nx - 2) + m], p_h_R), UL_y[5*(Ny - 1)*(Nx - 2) + m], 
			UR_y[5*(Ny - 1)*(Nx - 2) + m], CurY_starL_h, CurY_starR_h, SL_h, SR_h, SM_h);

		if (num_vars > 6){
			//Energy density of intermediate state
			float E_starL_e = ((SL_e - vy_e_L)*UL_y[6*(Ny - 1)*(Nx - 2) + m] + p_e_L*(SM_e - vy_e_L) + SM_e*UL_y[m]*(SM_e - vy_e_L)*(SL_e - vy_e_L))/(SL_e - SM_e);
			float E_starR_e = ((SR_e - vy_e_R)*UR_y[6*(Ny - 1)*(Nx - 2) + m] + p_e_R*(SM_e - vy_e_R) + SM_e*UR_y[m]*(SM_e - vy_e_R)*(SR_e - vy_e_R))/(SR_e - SM_e);
			float E_starL_h = ((SL_h - vy_h_L)*UL_y[7*(Ny - 1)*(Nx - 2) + m] + p_h_L*(SM_h - vy_h_L) + SM_h*UL_y[(Ny - 1)*(Nx - 2) + m]*(SM_h - vy_h_L)*(SL_h - vy_h_L))/(SL_h - SM_h);
			float E_starR_h = ((SR_h - vy_h_R)*UR_y[7*(Ny - 1)*(Nx - 2) + m] + p_h_R*(SM_h - vy_h_R) + SM_h*UR_y[(Ny - 1)*(Nx - 2) + m]*(SM_h - vy_h_R)*(SR_h - vy_h_R))/(SR_h - SM_h);

			//Electron internal energy
			G[6*(Ny - 1)*(Nx - 2) + m] = HLLFlux(EnergyFluxY2species(UL_y[m], UL_y[3*(Ny - 1)*(Nx - 2) + m], UL_y[6*(Ny - 1)*(Nx - 2) + m], p_e_L), 
				EnergyFluxY2species(UR_y[m], UR_y[3*(Ny - 1)*(Nx - 2) + m], UR_y[6*(Ny - 1)*(Nx - 2) + m], p_e_R), 
				UL_y[6*(Ny - 1)*(Nx - 2) + m], UR_y[6*(Ny - 1)*(Nx - 2) + m], E_starL_e, E_starR_e, SL_e, SR_e, SM_e);

			//Hole internal energy
			G[7*(Ny - 1)*(Nx - 2) + m] = HLLFlux(EnergyFluxY2species(UL_y[(Ny - 1)*(Nx - 2) + m], UL_y[5*(Ny - 1)*(Nx - 2) + m], UL_y[7*(Ny - 1)*(Nx - 2) + m], p_h_L), 
				EnergyFluxY2species(UR_y[(Ny - 1)*(Nx - 2) + m], UR_y[5*(Ny - 1)*(Nx - 2) + m], UR_y[7*(Ny - 1)*(Nx - 2) + m], p_h_R), 
				UL_y[7*(Ny - 1)*(Nx - 2) + m], UR_y[7*(Ny - 1)*(Nx - 2) + m], E_starL_h, E_starR_h, SL_h, SR_h, SM_h);
		}
	}	
}

void GrapheneCNP2D::MUSCL_RK3(){
	//Reconstruction applying a slope limiter for each of the variables
	Reconstruction(Den, Den_h, CurX, CurY, CurX_h, CurY_h, E, E_h, P, P_h, Tmp, Tmp_h);		
	
	//Compute the flux for the initial time step
	ComputeFluxX();
	ComputeFluxY();

	CflCondition();

	time_counter += dt;
	cout << vmax_x << "\t" << vmax_y << "\t" << time_counter << "/" << this->GetTmax() << endl;

	//
	//First time-step
	//

#pragma omp parallel for  default(none) shared(F, G, den_half, den_h_half)
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			den_half[m] = Den[m] - (dt/dx)*(F[m_mid_x] - F[m_mid_x - 1]) - (dt/dy)*(G[m_mid_y] - G[m_mid_y - 1]);
			
			den_h_half[m] = Den_h[m] - (dt/dx)*(F[(Nx - 1)*(Ny - 2) + m_mid_x] - F[(Nx - 1)*(Ny - 2) + m_mid_x - 1]) 
			- (dt/dy)*(G[(Ny - 1)*(Nx - 2) + m_mid_y] - G[(Ny - 1)*(Nx - 2) + m_mid_y - 1]);
		}
		den_half[i] = den_half[i + (Ny - 2)*Nx];
		den_half[i + (Ny - 1)*Nx] = den_half[i + Nx]; 
		
		den_h_half[i] = den_h_half[i + (Ny - 2)*Nx];
		den_h_half[i + (Ny - 1)*Nx] = den_h_half[i + Nx]; 
	}

	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		den_half[m0] = den_half[m0 + Nx - 2];
		den_half[mN] = den_half[m0 + 1];
		
		den_h_half[m0] = den_h_half[m0 + Nx - 2];
		den_h_half[mN] = den_h_half[m0 + 1];
	}

	//Compute the electric force at the first step and use an implicit source term
	GetPhi2species(den_half, den_h_half);

#pragma omp parallel for  default(none) shared(F, G, curX_half, curY_half, curX_h_half, curY_h_half, EX, EY, E0, e_half, e_h_half)
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			float F1 = CurX[m] - (dt/dx)*(F[2*(Nx - 1)*(Ny - 2) + m_mid_x] - F[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(dt/dy)*(G[2*(Ny - 1)*(Nx - 2) + m_mid_y] - G[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - dt*den_half[m]*(EX[m] + E0);

			float F2 = CurY[m] - (dt/dx)*(F[3*(Nx - 1)*(Ny - 2) + m_mid_x] - F[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(dt/dy)*(G[3*(Ny - 1)*(Nx - 2) + m_mid_y] - G[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - dt*den_half[m]*EY[m];

			float F3 = CurX_h[m] - (dt/dx)*(F[4*(Nx - 1)*(Ny - 2) + m_mid_x] - F[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(dt/dy)*(G[4*(Ny - 1)*(Nx - 2) + m_mid_y] - G[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + dt*den_h_half[m]*(EX[m] + E0);

			float F4 = CurY_h[m] - (dt/dx)*(F[5*(Nx - 1)*(Ny - 2) + m_mid_x] - F[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(dt/dy)*(G[5*(Ny - 1)*(Nx - 2) + m_mid_y] - G[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + dt*den_h_half[m]*EY[m];

			float r = den_half[m]/den_h_half[m];
			float r_e = den_h_half[m]/(den_half[m] + den_h_half[m]);
			float r_h = den_half[m]/(den_half[m] + den_h_half[m]);

			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*Tmp[m];
			}

			curX_half[m] = (F1 + dt*(-F2*cyc_freq + F3*r*3.*col_freq*r_e) + dt*F1*(3.*col_freq*r_e + col_freq))/(1. + 2.*dt*(3.*col_freq*r_e + col_freq) 
				+ dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq)));

			curY_half[m] = (F2 + dt*(F1*cyc_freq + F4*r*3.*col_freq*r_e) + dt*F2*(3.*col_freq*r_e + col_freq))/(1. + 2.*dt*(3.*col_freq*r_e + col_freq) 
				+ dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq))); 

			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*Tmp_h[m];
			}
			
			curX_h_half[m] = (F3 + dt*(F4*cyc_freq + F1*(1./r)*3.*col_freq*r_h) + dt*F3*(3.*col_freq*r_h + col_freq))/(1. + 2.*dt*(3.*col_freq*r_h + col_freq) 
				+ dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			curY_h_half[m] = (F4 + dt*(-F3*cyc_freq + F2*(1./r)*3.*col_freq*r_h) + dt*F4*(3.*col_freq*r_h + col_freq))/(1. + 2.*dt*(3.*col_freq*r_h + col_freq) 
				+ dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			float vel_x = curX_half[m]/den_half[m];
			float vel_y = curY_half[m]/den_half[m];
			float vel_h_x = curX_h_half[m]/den_h_half[m];
			float vel_h_y = curY_h_half[m]/den_h_half[m];

			if (num_vars > 6){
				if (col_freq > 0.){
					col_freq = 1.2522*Tmp[m];
				}
				e_half[m] = E[m] - (dt/dx)*(F[6*(Nx - 1)*(Ny - 2) + m_mid_x] - F[6*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				(dt/dy)*(G[6*(Ny - 1)*(Nx - 2) + m_mid_y] - G[6*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - dt*(curX_half[m]*(EX[m] + E0) + curY_half[m]*EY[m])
				- dt*(0.0697113/log(2.))*log(sqrt(exp(den_half[m]*log(2.)/Tmp[m]) - 1.) 
					+ 1./sqrt(exp(den_half[m]*log(2.)/Tmp[m]) - 1.))*col_freq*(Tmp[m] - save_T_e);
				//- dt*(curX_half[m]*(col_freq*vel_x + r_e*3.*col_freq*(vel_x - vel_h_x)) + curY_half[m]*(col_freq*vel_y + r_e*3.*col_freq*(vel_y - vel_h_y)));

				if (col_freq > 0.){
					col_freq = 1.2522*Tmp_h[m];
				}
				e_h_half[m] = E_h[m] - (dt/dx)*(F[7*(Nx - 1)*(Ny - 2) + m_mid_x] - F[7*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				(dt/dy)*(G[7*(Ny - 1)*(Nx - 2) + m_mid_y] - G[7*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + dt*(curX_h_half[m]*(EX[m] + E0) + curY_h_half[m]*EY[m])
				- dt*(0.0697113/log(2.))*log(sqrt(exp(den_h_half[m]*log(2.)/Tmp_h[m]) - 1.) 
					+ 1./sqrt(exp(den_h_half[m]*log(2.)/Tmp_h[m]) - 1.))*col_freq*(Tmp_h[m] - save_T_h);
				//- dt*(curX_h_half[m]*(col_freq*vel_h_x + r_h*3.*col_freq*(vel_h_x - vel_x)) + curY_h_half[m]*(col_freq*vel_h_y + r_h*3.*col_freq*(vel_h_y - vel_y)));

				p_half[m] = e_half[m] - 0.5*den_half[m]*(vel_x*vel_x + vel_y*vel_y);
				if (e_half[m] < 0.){
					printf("e = %f \n", e_half[m]);
				}
				p_h_half[m] = e_h_half[m] - 0.5*den_h_half[m]*(vel_h_x*vel_h_x + vel_h_y*vel_h_y);
				if (e_h_half[m] < 0.){
					printf("e = %f \n", e_h_half[m]);
				}
				t_half[m] = GetTmp(p_half[m], den_half[m], Tmp[m]);
				t_h_half[m] = GetTmp(p_h_half[m], den_h_half[m], Tmp_h[m]);
			}
		}
		//Boundary conditions along y
		curX_half[i] = curX_half[i + (Ny - 2)*Nx];
		curX_half[i + (Ny - 1)*Nx] = curX_half[i + Nx]; 
		
		curY_half[i] = curY_half[i + (Ny - 2)*Nx];
		curY_half[i + (Ny - 1)*Nx] = curY_half[i + Nx]; 

		curX_h_half[i] = curX_h_half[i + (Ny - 2)*Nx];
		curX_h_half[i + (Ny - 1)*Nx] = curX_h_half[i + Nx]; 
		
		curY_h_half[i] = curY_h_half[i + (Ny - 2)*Nx];
		curY_h_half[i + (Ny - 1)*Nx] = curY_h_half[i + Nx]; 

		if (num_vars > 6){
			e_half[i] = e_half[i + (Ny - 2)*Nx];
			e_half[i + (Ny - 1)*Nx] = e_half[i + Nx];

			e_h_half[i] = e_h_half[i + (Ny - 2)*Nx];
			e_h_half[i + (Ny - 1)*Nx] = e_h_half[i + Nx];

			p_half[i] = p_half[i + (Ny - 2)*Nx];
			p_half[i + (Ny - 1)*Nx] = p_half[i + Nx];

			p_h_half[i] = p_h_half[i + (Ny - 2)*Nx];
			p_h_half[i + (Ny - 1)*Nx] = p_h_half[i + Nx];

			t_half[i] = t_half[i + (Ny - 2)*Nx];
			t_half[i + (Ny - 1)*Nx] = t_half[i + Nx];

			t_h_half[i] = t_h_half[i + (Ny - 2)*Nx];
			t_h_half[i + (Ny - 1)*Nx] = t_h_half[i + Nx];
		}
	}

	//Boundary conditions along x
	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		curX_half[m0] = curX_half[m0 + Nx - 2];
		curX_half[mN] = curX_half[m0 + 1]; 
		
		curY_half[m0] = curY_half[m0 + Nx - 2];
		curY_half[mN] = curY_half[m0 + 1]; 

		curX_h_half[m0] = curX_h_half[m0 + Nx - 2];
		curX_h_half[mN] = curX_h_half[m0 + 1]; 
		
		curY_h_half[m0] = curY_h_half[m0 + Nx - 2];
		curY_h_half[mN] = curY_h_half[m0 + 1]; 

		if (num_vars > 6){
			e_half[m0] = e_half[m0 + Nx - 2];
			e_half[mN] = e_half[m0 + 1];

			e_h_half[m0] = e_h_half[m0 + Nx - 2];
			e_h_half[mN] = e_h_half[m0 + 1];

			p_half[m0] = p_half[m0 + Nx - 2];
			p_half[mN] = p_half[m0 + 1];

			p_h_half[m0] = p_h_half[m0 + Nx - 2];
			p_h_half[mN] = p_h_half[m0 + 1];

			t_half[m0] = t_half[m0 + Nx - 2];
			t_half[mN] = t_half[m0 + 1];

			t_h_half[m0] = t_h_half[m0 + Nx - 2];
			t_h_half[mN] = t_h_half[m0 + 1];
		}
	}

	//Repeat the procedure for the new values
	Reconstruction(den_half, den_h_half, curX_half, curY_half, curX_h_half, curY_h_half, e_half, e_h_half, p_half, p_h_half, t_half, t_h_half);

	ComputeFluxX();
	ComputeFluxY();

	//
	//Second time-step
	//

#pragma omp parallel for  default(none) shared(F, G, den_half, den_h_half)
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			den_half[m] = 0.75*Den[m] + 0.25*den_half[m] - (0.25*dt/dx)*(F[m_mid_x] - F[m_mid_x - 1]) - (0.25*dt/dy)*(G[m_mid_y] - G[m_mid_y - 1]);
			
			den_h_half[m] = 0.75*Den_h[m] + 0.25*den_h_half[m] - (0.25*dt/dx)*(F[(Nx - 1)*(Ny - 2) + m_mid_x] - F[(Nx - 1)*(Ny - 2) + m_mid_x - 1]) 
			- (0.25*dt/dy)*(G[(Ny - 1)*(Nx - 2) + m_mid_y] - G[(Ny - 1)*(Nx - 2) + m_mid_y - 1]);
		}
		den_half[i] = den_half[i + (Ny - 2)*Nx];
		den_half[i + (Ny - 1)*Nx] = den_half[i + Nx]; 
		
		den_h_half[i] = den_h_half[i + (Ny - 2)*Nx];
		den_h_half[i + (Ny - 1)*Nx] = den_h_half[i + Nx]; 
	}

	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		den_half[m0] = den_half[m0 + Nx - 2];
		den_half[mN] = den_half[m0 + 1];
		den_h_half[m0] = den_h_half[m0 + Nx - 2];
		den_h_half[mN] = den_h_half[m0 + 1];
	}

	//Compute the potential at the full step and use an implicit source term
	GetPhi2species(den_half, den_h_half);

#pragma omp parallel for  default(none) shared(F, G, curX_half, curY_half, curX_h_half, curY_h_half, EX, EY, E0, e_half, e_h_half)
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			float F1 = 0.75*CurX[m] + 0.25*curX_half[m] - (0.25*dt/dx)*(F[2*(Nx - 1)*(Ny - 2) + m_mid_x] - F[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(0.25*dt/dy)*(G[2*(Ny - 1)*(Nx - 2) + m_mid_y] - G[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - 0.25*dt*den_half[m]*(EX[m] + E0);

			float F2 = 0.75*CurY[m] + 0.25*curY_half[m] - (0.25*dt/dx)*(F[3*(Nx - 1)*(Ny - 2) + m_mid_x] - F[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(0.25*dt/dy)*(G[3*(Ny - 1)*(Nx - 2) + m_mid_y] - G[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - 0.25*dt*den_half[m]*EY[m];

			float F3 = 0.75*CurX_h[m] + 0.25*curX_h_half[m] - (0.25*dt/dx)*(F[4*(Nx - 1)*(Ny - 2) + m_mid_x] - F[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(0.25*dt/dy)*(G[4*(Ny - 1)*(Nx - 2) + m_mid_y] - G[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + 0.25*dt*den_h_half[m]*(EX[m] + E0);

			float F4 = 0.75*CurY_h[m] + 0.25*curY_h_half[m] - (0.25*dt/dx)*(F[5*(Nx - 1)*(Ny - 2) + m_mid_x] - F[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			(0.25*dt/dy)*(G[5*(Ny - 1)*(Nx - 2) + m_mid_y] - G[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + 0.25*dt*den_h_half[m]*EY[m];

			float r = den_half[m]/den_h_half[m];
			float r_e = den_h_half[m]/(den_half[m] + den_h_half[m]);
			float r_h = den_half[m]/(den_half[m] + den_h_half[m]);

			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*t_half[m];
			}
			
			curX_half[m] = (F1 + 0.25*dt*(-F2*cyc_freq + F3*r*3.*col_freq*r_e) + 0.25*dt*F1*(3.*col_freq*r_e + col_freq))/(1. + 2.*0.25*dt*(3.*col_freq*r_e + col_freq) 
				+ 0.25*0.25*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq)));

			curY_half[m] = (F2 + 0.25*dt*(F1*cyc_freq + F4*r*3.*col_freq*r_e) + 0.25*dt*F2*(3.*col_freq*r_e + col_freq))/(1. + 2.*0.25*dt*(3.*col_freq*r_e + col_freq) 
				+ 0.25*0.25*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq))); 

			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*t_h_half[m];
			}
			curX_h_half[m] = (F3 + 0.25*dt*(F4*cyc_freq + F1*(1./r)*3.*col_freq*r_h) + 0.25*dt*F3*(3.*col_freq*r_h + col_freq))/(1. + 2.*0.25*dt*(3.*col_freq*r_h + col_freq) 
				+ 0.25*0.25*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			curY_h_half[m] = (F4 + 0.25*dt*(-F3*cyc_freq + F2*(1./r)*3.*col_freq*r_h) + 0.25*dt*F4*(3.*col_freq*r_h + col_freq))/(1. + 2.*0.25*dt*(3.*col_freq*r_h + col_freq) 
				+ 0.25*0.25*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			float vel_x = curX_half[m]/den_half[m];
			float vel_y = curY_half[m]/den_half[m];
			float vel_h_x = curX_h_half[m]/den_h_half[m];
			float vel_h_y = curY_h_half[m]/den_h_half[m];

			if (num_vars > 6){
				if (col_freq > 0.){
					col_freq = 1.2522*t_half[m];
				}
				e_half[m] = 0.75*E[m] + 0.25*e_half[m] - (0.25*dt/dx)*(F[6*(Nx - 1)*(Ny - 2) + m_mid_x] - F[6*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				(0.25*dt/dy)*(G[6*(Ny - 1)*(Nx - 2) + m_mid_y] - G[6*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - 0.25*dt*(curX_half[m]*(EX[m] + E0) + curY_half[m]*EY[m])
				- 0.25*dt*(0.0697113/log(2.))*log(sqrt(exp(den_half[m]*log(2.)/t_half[m]) - 1.) 
					+ 1./sqrt(exp(den_half[m]*log(2.)/t_half[m]) - 1.))*col_freq*(t_half[m] - save_T_e); 
				//- 0.25*dt*(curX_half[m]*(col_freq*vel_x + r_e*3.*col_freq*(vel_x - vel_h_x)) + curY_half[m]*(col_freq*vel_y + r_e*3.*col_freq*(vel_y - vel_h_y)));

				if (col_freq > 0.){
					col_freq = 1.2522*t_h_half[m];
				}
				e_h_half[m] = 0.75*E_h[m] + 0.25*e_h_half[m] - (0.25*dt/dx)*(F[7*(Nx - 1)*(Ny - 2) + m_mid_x] - F[7*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				(0.25*dt/dy)*(G[7*(Ny - 1)*(Nx - 2) + m_mid_y] - G[7*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + 0.25*dt*(curX_h_half[m]*(EX[m] + E0) + curY_h_half[m]*EY[m])
				- 0.25*dt*(0.0697113/log(2.))*log(sqrt(exp(den_h_half[m]*log(2.)/t_h_half[m]) - 1.) 
					+ 1./sqrt(exp(den_h_half[m]*log(2.)/t_h_half[m]) - 1.))*col_freq*(t_h_half[m] - save_T_h);
				//- 0.25*dt*(curX_h_half[m]*(col_freq*vel_h_x + r_h*3.*col_freq*(vel_h_x - vel_x)) + curY_h_half[m]*(col_freq*vel_h_y + r_h*3.*col_freq*(vel_h_y - vel_y)));

				p_half[m] = e_half[m] - 0.5*den_half[m]*(vel_x*vel_x + vel_y*vel_y);
				if (e_half[m] < 0.){
					printf("e = %f \n", e_half[m]);
				}
				p_h_half[m] = e_h_half[m] - 0.5*den_h_half[m]*(vel_h_x*vel_h_x + vel_h_y*vel_h_y);
				if (e_h_half[m] < 0.){
					printf("e = %f \n", e_h_half[m]);
				}
				t_half[m] = GetTmp(p_half[m], den_half[m], t_half[m]);
				t_h_half[m] = GetTmp(p_h_half[m], den_h_half[m], t_h_half[m]);
			}
		}
		//Boundary conditions along y
		curX_half[i] = curX_half[i + (Ny - 2)*Nx];
		curX_half[i + (Ny - 1)*Nx] = curX_half[i + Nx]; 
		
		curY_half[i] = curY_half[i + (Ny - 2)*Nx];
		curY_half[i + (Ny - 1)*Nx] = curY_half[i + Nx]; 

		curX_h_half[i] = curX_h_half[i + (Ny - 2)*Nx];
		curX_h_half[i + (Ny - 1)*Nx] = curX_h_half[i + Nx]; 
		
		curY_h_half[i] = curY_h_half[i + (Ny - 2)*Nx];
		curY_h_half[i + (Ny - 1)*Nx] = curY_h_half[i + Nx];

		if (num_vars > 6){
			e_half[i] = e_half[i + (Ny - 2)*Nx];
			e_half[i + (Ny - 1)*Nx] = e_half[i + Nx];

			e_h_half[i] = e_h_half[i + (Ny - 2)*Nx];
			e_h_half[i + (Ny - 1)*Nx] = e_h_half[i + Nx];

			p_half[i] = p_half[i + (Ny - 2)*Nx];
			p_half[i + (Ny - 1)*Nx] = p_half[i + Nx];

			p_h_half[i] = p_h_half[i + (Ny - 2)*Nx];
			p_h_half[i + (Ny - 1)*Nx] = p_h_half[i + Nx];

			t_half[i] = t_half[i + (Ny - 2)*Nx];
			t_half[i + (Ny - 1)*Nx] = t_half[i + Nx];

			t_h_half[i] = t_h_half[i + (Ny - 2)*Nx];
			t_h_half[i + (Ny - 1)*Nx] = t_h_half[i + Nx];
		} 
	}

	//Boundary conditions along x
	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		curX_half[m0] = curX_half[m0 + Nx - 2];
		curX_half[mN] = curX_half[m0 + 1]; 
		
		curY_half[m0] = curY_half[m0 + Nx - 2];
		curY_half[mN] = curY_half[m0 + 1]; 

		curX_h_half[m0] = curX_h_half[m0 + Nx - 2];
		curX_h_half[mN] = curX_h_half[m0 + 1]; 
		
		curY_h_half[m0] = curY_h_half[m0 + Nx - 2];
		curY_h_half[mN] = curY_h_half[m0 + 1]; 

		if (num_vars > 6){
			e_half[m0] = e_half[m0 + Nx - 2];
			e_half[mN] = e_half[m0 + 1];

			e_h_half[m0] = e_h_half[m0 + Nx - 2];
			e_h_half[mN] = e_h_half[m0 + 1];

			p_half[m0] = p_half[m0 + Nx - 2];
			p_half[mN] = p_half[m0 + 1];

			p_h_half[m0] = p_h_half[m0 + Nx - 2];
			p_h_half[mN] = p_h_half[m0 + 1];

			t_half[m0] = t_half[m0 + Nx - 2];
			t_half[mN] = t_half[m0 + 1];

			t_h_half[m0] = t_h_half[m0 + Nx - 2];
			t_h_half[mN] = t_h_half[m0 + 1];
		}
	}

	//Repeat the procedure for the new values
	Reconstruction(den_half, den_h_half, curX_half, curY_half, curX_h_half, curY_h_half, e_half, e_h_half, p_half, p_h_half, t_half, t_h_half);

	ComputeFluxX();
	ComputeFluxY();

	//
	//Third time-step
	//
#pragma omp parallel for  default(none) shared(F, G, Den, Den_h)
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			Den[m] = (1./3.)*Den[m] + (2./3.)*den_half[m] - ((2./3.)*dt/dx)*(F[m_mid_x] - F[m_mid_x - 1]) - ((2./3.)*dt/dy)*(G[m_mid_y] - G[m_mid_y - 1]);
			
			Den_h[m] = (1./3.)*Den_h[m] + (2./3.)*den_h_half[m] - ((2./3.)*dt/dx)*(F[(Nx - 1)*(Ny - 2) + m_mid_x] - F[(Nx - 1)*(Ny - 2) + m_mid_x - 1]) 
			- ((2./3.)*dt/dy)*(G[(Ny - 1)*(Nx - 2) + m_mid_y] - G[(Ny - 1)*(Nx - 2) + m_mid_y - 1]);
		}
		Den[i] = Den[i + (Ny - 2)*Nx];
		Den[i + (Ny - 1)*Nx] = Den[i + Nx]; 
		
		Den_h[i] = Den_h[i + (Ny - 2)*Nx];
		Den_h[i + (Ny - 1)*Nx] = Den_h[i + Nx]; 
	}

	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		Den[m0] = Den[m0 + Nx - 2];
		Den[mN] = Den[m0 + 1];
		
		Den_h[m0] = Den_h[m0 + Nx - 2];
		Den_h[mN] = Den_h[m0 + 1];
	}

	//Compute the potential at the full step and use an implicit source term
	GetPhi2species(Den, Den_h);

#pragma omp parallel for  default(none) shared(F, G, CurX, CurY, CurX_h, CurY_h, EX, EY, E0, E, E_h)
	//Add implicit magnetic field and collision source terms
	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;
			int m_mid_x = i + (j - 1)*(Nx - 1);
			int m_mid_y = j + (i - 1)*(Ny - 1);

			float F1 = (1./3.)*CurX[m] + (2./3.)*curX_half[m] - ((2./3.)*dt/dx)*(F[2*(Nx - 1)*(Ny - 2) + m_mid_x] - F[2*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			((2./3.)*dt/dy)*(G[2*(Ny - 1)*(Nx - 2) + m_mid_y] - G[2*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - (2./3.)*dt*Den[m]*(EX[m] + E0);

			float F2 = (1./3.)*CurY[m] + (2./3.)*curY_half[m] - ((2./3.)*dt/dx)*(F[3*(Nx - 1)*(Ny - 2) + m_mid_x] - F[3*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			((2./3.)*dt/dy)*(G[3*(Ny - 1)*(Nx - 2) + m_mid_y] - G[3*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - (2./3.)*dt*Den[m]*EY[m];

			float F3 = (1./3.)*CurX_h[m] + (2./3.)*curX_h_half[m] - ((2./3.)*dt/dx)*(F[4*(Nx - 1)*(Ny - 2) + m_mid_x] - F[4*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			((2./3.)*dt/dy)*(G[4*(Ny - 1)*(Nx - 2) + m_mid_y] - G[4*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + (2./3.)*dt*Den_h[m]*(EX[m] + E0);

			float F4 = (1./3.)*CurY_h[m] + (2./3.)*curY_h_half[m] - ((2./3.)*dt/dx)*(F[5*(Nx - 1)*(Ny - 2) + m_mid_x] - F[5*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
			((2./3.)*dt/dy)*(G[5*(Ny - 1)*(Nx - 2) + m_mid_y] - G[5*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + (2./3.)*dt*Den_h[m]*EY[m];

			float r = Den[m]/Den_h[m];
			float r_e = Den_h[m]/(Den[m] + Den_h[m]);
			float r_h = Den[m]/(Den[m] + Den_h[m]);

			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*t_half[m];
			}
			CurX[m] = (F1 + (2./3.)*dt*(-F2*cyc_freq + F3*r*3.*col_freq*r_e) + (2./3.)*dt*F1*(3.*col_freq*r_e + col_freq))/(1. + 2.*(2./3.)*dt*(3.*col_freq*r_e + col_freq) 
				+ (2./3.)*(2./3.)*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq)));

			CurY[m] = (F2 + (2./3.)*dt*(F1*cyc_freq + F4*r*3.*col_freq*r_e) + (2./3.)*dt*F2*(3.*col_freq*r_e + col_freq))/(1. + 2.*(2./3.)*dt*(3.*col_freq*r_e + col_freq) 
				+ (2./3.)*(2./3.)*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_e + col_freq))); 


			if (num_vars > 6 && col_freq > 0.){
				col_freq = 1.2522*t_h_half[m];
			}

			CurX_h[m] = (F3 + (2./3.)*dt*(F4*cyc_freq + F1*(1./r)*3.*col_freq*r_h) + (2./3.)*dt*F3*(3.*col_freq*r_h + col_freq))/(1. + 2.*(2./3.)*dt*(3.*col_freq*r_h + col_freq) 
				+ (2./3.)*(2./3.)*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			CurY_h[m] = (F4 + (2./3.)*dt*(-F3*cyc_freq + F2*(1./r)*3.*col_freq*r_h) + (2./3.)*dt*F4*(3.*col_freq*r_h + col_freq))/(1. + 2.*(2./3.)*dt*(3.*col_freq*r_h + col_freq) 
				+ (2./3.)*(2./3.)*dt*dt*(cyc_freq*cyc_freq + col_freq*(3.*col_freq*r_h + col_freq)));

			float vel_x = CurX[m]/Den[m];
			float vel_y = CurY[m]/Den[m];
			float vel_h_x = CurX_h[m]/Den_h[m];
			float vel_h_y = CurY_h[m]/Den_h[m];

			if (num_vars > 6){
				if (col_freq > 0.){
					col_freq = 1.2522*t_half[m];
				}
				E[m] = (1./3.)*E[m] + (2./3.)*e_half[m] - ((2./3.)*dt/dx)*(F[6*(Nx - 1)*(Ny - 2) + m_mid_x] - F[6*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				((2./3.)*dt/dy)*(G[6*(Ny - 1)*(Nx - 2) + m_mid_y] - G[6*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) - (2./3.)*dt*(CurX[m]*(EX[m] + E0) + CurY[m]*EY[m])
				- (2./3.)*dt*(0.0697113/log(2.))*log(sqrt(exp(Den[m]*log(2.)/t_half[m]) - 1.) 
					+ 1./sqrt(exp(Den[m]*log(2.)/t_half[m]) - 1.))*col_freq*(t_half[m] - save_T_e);
				//- (2./3.)*dt*(CurX[m]*(col_freq*vel_x + r_e*3.*col_freq*(vel_x - vel_h_x)) + CurY[m]*(col_freq*vel_y + r_e*3.*col_freq*(vel_y - vel_h_y)));

				if (col_freq > 0.){
					col_freq = 1.2522*t_h_half[m];
				}
				E_h[m] = (1./3.)*E_h[m] + (2./3.)*e_h_half[m] - ((2./3.)*dt/dx)*(F[7*(Nx - 1)*(Ny - 2) + m_mid_x] - F[7*(Nx - 1)*(Ny - 2) + m_mid_x - 1]) -
				((2./3.)*dt/dy)*(G[7*(Ny - 1)*(Nx - 2) + m_mid_y] - G[7*(Ny - 1)*(Nx - 2) + m_mid_y - 1]) + (2./3.)*dt*(CurX_h[m]*(EX[m] + E0) + CurY_h[m]*EY[m])
				- (2./3.)*dt*(0.0697113/log(2.))*log(sqrt(exp(Den_h[m]*log(2.)/t_h_half[m]) - 1.) 
					+ 1./sqrt(exp(Den_h[m]*log(2.)/t_h_half[m]) - 1.))*col_freq*(t_h_half[m] - save_T_h);
				//- (2./3.)*dt*(CurX_h[m]*(col_freq*vel_h_x + r_h*3.*col_freq*(vel_h_x - vel_x)) + CurY_h[m]*(col_freq*vel_h_y + r_h*3.*col_freq*(vel_h_y - vel_y)));

				P[m] = E[m] - 0.5*Den[m]*(vel_x*vel_x + vel_y*vel_y);
				if (E[m] < 0.){
					printf("e = %f \n", E[m]);
				}
				P_h[m] = E_h[m] - 0.5*Den_h[m]*(vel_h_x*vel_h_x + vel_h_y*vel_h_y);
				if (E_h[m] < 0.){
					printf("e = %f \n", E_h[m]);
				}
				Tmp[m] = GetTmp(P[m], Den[m], 1.0f);
				Tmp_h[m] = GetTmp(P_h[m], Den_h[m], 1.0f);
			}
		}
		CurX[i] = CurX[i + (Ny - 2)*Nx];
		CurX[i + (Ny - 1)*Nx] = CurX[i + Nx]; 
		
		CurY[i] = CurY[i + (Ny - 2)*Nx];
		CurY[i + (Ny - 1)*Nx] = CurY[i + Nx]; 

		CurX_h[i] = CurX_h[i + (Ny - 2)*Nx];
		CurX_h[i + (Ny - 1)*Nx] = CurX_h[i + Nx]; 
		
		CurY_h[i] = CurY_h[i + (Ny - 2)*Nx];
		CurY_h[i + (Ny - 1)*Nx] = CurY_h[i + Nx]; 

		if (num_vars > 6){
			E[i] = E[i + (Ny - 2)*Nx];
			E[i + (Ny - 1)*Nx] = E[i + Nx];

			E_h[i] = E_h[i + (Ny - 2)*Nx];
			E_h[i + (Ny - 1)*Nx] = E_h[i + Nx];

			P[i] = P[i + (Ny - 2)*Nx];
			P[i + (Ny - 1)*Nx] = P[i + Nx];

			P_h[i] = P_h[i + (Ny - 2)*Nx];
			P_h[i + (Ny - 1)*Nx] = P_h[i + Nx];

			Tmp[i] = Tmp[i + (Ny - 2)*Nx];
			Tmp[i + (Ny - 1)*Nx] = Tmp[i + Nx];

			Tmp_h[i] = Tmp_h[i + (Ny - 2)*Nx];
			Tmp_h[i + (Ny - 1)*Nx] = Tmp_h[i + Nx];
		}

		EX[i] = EX[i + (Ny - 2)*Nx];
		EX[i + (Ny - 1)*Nx] = EX[i + Nx];
		EY[i] = EY[i + (Ny - 2)*Nx];
		EY[i + (Ny - 1)*Nx] = EY[i + Nx];
	}

	//Boundary conditions along x
	for (int j = 0; j < Ny; ++j){
		int m0 = j*Nx;
		int mN = Nx - 1 + j*Nx;
		CurX[m0] = CurX[m0 + Nx - 2];
		CurX[mN] = CurX[m0 + 1]; 
		
		CurY[m0] = CurY[m0 + Nx - 2];
		CurY[mN] = CurY[m0 + 1]; 

		CurX_h[m0] = CurX_h[m0 + Nx - 2];
		CurX_h[mN] = CurX_h[m0 + 1]; 
		
		CurY_h[m0] = CurY_h[m0 + Nx - 2];
		CurY_h[mN] = CurY_h[m0 + 1];

		if (num_vars > 6){
			E[m0] = E[m0 + Nx - 2];
			E[mN] = E[m0 + 1];

			E_h[m0] = E_h[m0 + Nx - 2];
			E_h[mN] = E_h[m0 + 1];

			P[m0] = P[m0 + Nx - 2];
			P[mN] = P[m0 + 1];

			P_h[m0] = P_h[m0 + Nx - 2];
			P_h[mN] = P_h[m0 + 1];

			Tmp[m0] = Tmp[m0 + Nx - 2];
			Tmp[mN] = Tmp[m0 + 1];

			Tmp_h[m0] = Tmp_h[m0 + Nx - 2];
			Tmp_h[mN] = Tmp_h[m0 + 1];
		}

		EX[m0] = EX[m0 + Nx - 2];
		EX[mN] = EX[m0 + 1]; 
		
		EY[m0] = EY[m0 + Nx - 2];
		EY[mN] = EY[m0 + 1]; 
	}
	
	vmax_x = 0.;
	vmax_y = 0.;

	float aux_e = 0.;
	float aux_h = 0.;
	float a = 0.;
	float b = 0.;

	/*float Cur_avgX = 0.;
	float Cur_avgY = 0.;
	for (int i = 0; i < Nx*Ny; ++i){
		Cur_avgX += CurX[i]/(Nx*Ny);
		Cur_avgY += CurY[i]/(Nx*Ny);
	}*/

	for (int i = 0; i < Nx*Ny; ++i){
		VelX[i] = CurX[i]/Den[i];
		VelY[i] = CurY[i]/Den[i];
		VelX_h[i] = CurX_h[i]/Den_h[i];
		VelY_h[i] = CurY_h[i]/Den_h[i];

		if (aux_e < Den[i]){
			aux_e = Den[i];
		}

		if (aux_h < Den_h[i]){
			aux_h = Den_h[i];
		}

		if (fabs(EX[i]) > a){
			a = fabs(EX[i]);
		}

		if (fabs(EY[i]) > b){
			b = fabs(EY[i]);
		}
		
		float cs_e = sqrt(log(2.)*Den[i]/(1. - exp(-log(2.)*Den[i])));
		float cs_h = sqrt(log(2.)*Den_h[i]/(1. - exp(-log(2.)*Den_h[i])));

		if (num_vars > 6){
			cs_e = sqrt(log(2.)*Den[i]/(1. - exp(-log(2.)*Den[i]/Tmp[i])));
			cs_h = sqrt(log(2.)*Den_h[i]/(1. - exp(-log(2.)*Den_h[i]/Tmp_h[i])));
		}

		float eig_max_x = max(fabs(VelX[i]) + cs_e, fabs(VelX_h[i]) + cs_h);
		float eig_max_y = max(fabs(VelY[i]) + cs_e, fabs(VelY_h[i]) + cs_h);

		if (vmax_x < eig_max_x){
			vmax_x = eig_max_x;
		}

		if (vmax_y < eig_max_y){
			vmax_y = eig_max_y;
		}

		den_E_kin[i] = 0.5*(Den[i]*(VelX[i]*VelX[i] + VelY[i]*VelY[i]) + Den_h[i]*(VelX_h[i]*VelX_h[i] + VelY_h[i]*VelY_h[i]));
		den_E_th[i] = -(save_T_e*save_T_e/log(2.))*Li2(1. - exp(log(2.)*Den[i]/save_T_e)) 
		- (save_T_h*save_T_h/log(2.))*Li2(1. - exp(log(2.)*Den_h[i]/save_T_h));

		if (num_vars > 6){
			den_E_th[i] = E[i] + E_h[i] - den_E_kin[i];
		}

		//den_E_energy[i] = 0.5*(EX[i]*EX[i] + EY[i]*EY[i]);
		cur_totX[i] = CurX_h[i] - CurX[i];
		cur_totY[i] = CurY_h[i] - CurY[i];
		den_E_energy[i] = EX[i]*CurX[i] + EY[i]*CurY[i];
	}
	J_power = Integral_2_D(Nx, Ny, dx, dy, den_E_energy);
	E_kin = Integral_2_D(Nx, Ny, dx, dy, den_E_kin);
	E_th = Integral_2_D(Nx, Ny, dx, dy, den_E_th);
	I_x = Integral_2_D(Nx, Ny, dx, dy, cur_totX);
	I_y = Integral_2_D(Nx, Ny, dx, dy, cur_totY);

	cout << aux_e << "\t" << E_th << "\t" << E_energy << "\t" << I_x << "\t" << I_y << endl << endl;
}

float GrapheneCNP2D::DensityFluxX2species(float cur){
	float f = cur;
	return f;
}

float GrapheneCNP2D::DensityFluxY2species(float cur){
	float f = cur;
	return f;
}

float GrapheneCNP2D::XCurrentFluxX2species(float den, float cur, float p){
	float f = cur*cur/den + p;
	return f;
}

float GrapheneCNP2D::XCurrentFluxY2species(float den, float cur1, float cur2){
	float f = cur1*cur2/den;
	return f;
}

float GrapheneCNP2D::YCurrentFluxX2species(float den, float cur1, float cur2){
	float f = cur1*cur2/den;
	return f;
}

float GrapheneCNP2D::YCurrentFluxY2species(float den, float cur, float p){
	float f = cur*cur/den + p;
	return f;
}

float GrapheneCNP2D::EnergyFluxX2species(float den, float cur, float e, float p){
	float f = (e + p)*cur/den;
	return f;
}

float GrapheneCNP2D::EnergyFluxY2species(float den, float cur, float e, float p){
	float f = (e + p)*cur/den;
	return f;
}

void GrapheneCNP2D::GetPhi2species(const float *n_e, const float *n_h){
	int n_size = Nx*Ny;
		
	fftwf_complex *n_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size);
	fftwf_complex *n_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size);
	fftwf_plan p_n = fftwf_plan_dft_2d(Nx, Ny, n_in, n_out, FFTW_FORWARD, FFTW_MEASURE);

	for (int i = 0; i < n_size; ++i){
		n_in[i][0] = n_h[i] - n_e[i];
		n_in[i][1] = 0.; 
	}
	
	fftwf_execute(p_n);

	fftwf_complex *V_q_x = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size); 
	fftwf_complex *V_r_x = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size); 
	fftwf_complex *V_q_y = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size); 
	fftwf_complex *V_r_y = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*n_size); 

	fftwf_plan p_V_x = fftwf_plan_dft_2d(Nx, Ny, V_q_x, V_r_x, FFTW_BACKWARD, FFTW_MEASURE);
	fftwf_plan p_V_y = fftwf_plan_dft_2d(Nx, Ny, V_q_y, V_r_y, FFTW_BACKWARD, FFTW_MEASURE);

	float constant = 1000.;
	float dqx = 2.*M_PI/(Nx*dx);
	float dqy = 2.*M_PI/(Ny*dy);

	E_energy = 0.;
	for (int i = 0; i < Nx; ++i){
		for (int j = 0; j < Ny; ++j){
			int m = i + j*Nx;
			if (i == 0 && j == 0){
				V_q_x[0][0] = 0.f;
				V_q_x[0][1] = 0.f;
				V_q_y[0][0] = 0.f;
				V_q_y[0][1] = 0.f;
				continue;
			}

			float qx;
			if (i <= Nx/2){
				qx = i*dqx;
			}else{
				qx = -(Nx - i)*dqx;
			}

			float qy;
			if (j <= Ny/2){
				qy = j*dqy;
			}
			else{
				qy = -(Ny - j)*dqy;
			}

			float norm = sqrt(qx*qx + qy*qy);

			if (norm <= dqx*64.f){
				V_q_x[m][1] = qx*(constant/norm)*n_out[m][0];
				V_q_x[m][0] = -qx*(constant/norm)*n_out[m][1];

				V_q_y[m][1] = qy*(constant/norm)*n_out[m][0];
				V_q_y[m][0] = -qy*(constant/norm)*n_out[m][1];

				E_energy += (constant*constant/norm)*(n_out[m][0]*n_out[m][0] + n_out[m][1]*n_out[m][1])/n_size;
			}
			else{
				V_q_x[m][1] = 0.f;
				V_q_x[m][0] = 0.f;
			
				V_q_y[m][1] = 0.f;
				V_q_y[m][0] = 0.f;
			}
		}
	}

	E_energy /= n_size;

	//aux_var = sqrt(V_q_x[56 + 30*Nx][1]*V_q_x[56 + 30*Nx][1] + V_q_x[56 + 30*Nx][0]*V_q_x[56 + 30*Nx][0])/n_size;

	fftwf_execute(p_V_x);
	fftwf_execute(p_V_y);

	for (int i = 1; i < Nx - 1; ++i){
		for (int j = 1; j < Ny - 1; ++j){
			int m = i + j*Nx;

			EX[m] = -V_r_x[m][0]/n_size;
			EY[m] = -V_r_y[m][0]/n_size;
		} 
	}

	fftwf_free(n_in);
	fftwf_free(n_out);
	fftwf_free(V_r_x);
	fftwf_free(V_q_x);
	fftwf_free(V_r_y);
	fftwf_free(V_q_y);
	fftwf_destroy_plan(p_n);
	fftwf_destroy_plan(p_V_x);
	fftwf_destroy_plan(p_V_y);
}