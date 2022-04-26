/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/DiracGraphene2DLib.h"



DiracGraphene2D::DiracGraphene2D(SetUpParameters &input_parameters) : Fluid2D(input_parameters) {
	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency;
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency;
	therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity
	char buffer [100];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
	file_infix = buffer;
}

void DiracGraphene2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

/*
void DiracGraphene2D::MassFluxToVelocity(const string grid) {
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

void DiracGraphene2D::CflCondition(){ // Eventual redefinition
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



float DiracGraphene2D::DensityFluxX(GridPoint2D p, char side) {
	float den;
	float px;

	px = SideAverage(ptr_px, p, side);
	den = SideAverage(ptr_den, p, side);

	return px / sqrt(den);  // DUVIDA: porque é que aqui nao se poe sobre a massa - DensityToMass?
}

float DiracGraphene2D::DensityFluxY(GridPoint2D p, char side) {
	float den;
	float py;
	py = SideAverage(ptr_py, p, side);
	den = SideAverage(ptr_den, p, side);
	return py / sqrt(den);  // DUVIDA: porque é que aqui nao se poe sobre a massa - DensityToMass?
}

float DiracGraphene2D::XMomentumFluxX(GridPoint2D p, char side) {

	float sound;
	float den;
	float hden;
	float px;
	float mass;

	sound = SideAverage(ptr_snd, p, side);
	den = SideAverage(ptr_den, p, side);
	hden = SideAverage(hptr_den, p, side);
	px = SideAverage(ptr_px, p, side);

	mass=DensityToMass(den);

	return px * px / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * (den - hden) ;
}

float DiracGraphene2D::XMomentumFluxY(GridPoint2D p, char side) {

	float den;
	float py;
	float px;
	float mass;



	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	mass=DensityToMass(den);

	return px * py / mass;
}


float DiracGraphene2D::YMomentumFluxY(GridPoint2D p, char side) {

	float sound ;
	float den;
	float hden;
	float py;
	float mass;

	sound = SideAverage(ptr_snd, p, side);
	den = SideAverage(ptr_den, p, side);
	hden = SideAverage(hptr_den, p, side);
	py = SideAverage(ptr_py, p, side);
	mass=DensityToMass(den);

	return py * py / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * (den - hden);
}

float DiracGraphene2D::YMomentumFluxX(GridPoint2D p, char side) {

	float den;
	float px;
	float py;
	float mass;

	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	mass=DensityToMass(den);

	return px * py / mass;
}


float DiracGraphene2D::DensitySource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

float DiracGraphene2D::XMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_x  - cyc_freq*flx_y/sqrt(n);
}

float DiracGraphene2D::YMomentumSource(float n, float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return -1.0f*col_freq*flx_y  + cyc_freq*flx_x/sqrt(n);
}


//________________________________ HOLE ____________________________

float DiracGraphene2D::HDensityFluxX(GridPoint2D p, char side) {
	float den;
	float px;

	px = SideAverage(hptr_px, p, side);
	den = SideAverage(hptr_den, p, side);

	return px / sqrt(den);
}

float DiracGraphene2D::HDensityFluxY(GridPoint2D p, char side) {
	float den;
	float py;
	py = SideAverage(hptr_py, p, side);
	den = SideAverage(hptr_den, p, side);
	return py / sqrt(den);
}

float DiracGraphene2D::HXMomentumFluxX(GridPoint2D p, char side) {

	float sound;
	float den;
	float hden;
	float px;
	float mass;

	sound = SideAverage(ptr_snd, p, side);
	hden = SideAverage(hptr_den, p, side);
	den = SideAverage(ptr_den, p, side);
	px = SideAverage(hptr_px, p, side);

	mass=DensityToMass(den);

	return px * px / mass + vel_fer * vel_fer * mass / 3.0f - 0.5f * sound * sound * (den - hden) ;
}

float DiracGraphene2D::HXMomentumFluxY(GridPoint2D p, char side) {

	float den;
	float py;
	float px;
	float mass;



	den = SideAverage(hptr_den, p, side);
	px = SideAverage(hptr_px, p, side);
	py = SideAverage(hptr_py, p, side);
	mass=DensityToMass(den);

	return px * py / mass;
}


float DiracGraphene2D::HYMomentumFluxY(GridPoint2D p, char side) {

	float sound ;
	float den;
	float hden;
	float py;
	float mass;

	sound = SideAverage(ptr_snd, p, side);
	hden = SideAverage(hptr_den, p, side);
	den = SideAverage(ptr_den, p, side);
	py = SideAverage(hptr_py, p, side);
	mass=DensityToMass(den);

	return py * py / mass + vel_fer * vel_fer * mass / 3.0f - 0.5f * sound * sound * (den - hden);
}

float DiracGraphene2D::HYMomentumFluxX(GridPoint2D p, char side) {

	float den;
	float px;
	float py;
	float mass;

	den = SideAverage(hptr_den, p, side);
	px = SideAverage(hptr_px, p, side);
	py = SideAverage(hptr_py, p, side);
	mass=DensityToMass(den);

	return px * py / mass;
}


DiracGraphene2D::~DiracGraphene2D(){
delete[] Den;
delete[] VelX;
delete[] VelY;
delete[] FlxX;
delete[] FlxY;
delete[] CurX;
delete[] CurY;
delete[] HDen;
delete[] HVelX;
delete[] HVelY;
delete[] HFlxX;
delete[] HFlxY;
delete[] HCurX;
delete[] HCurY;
delete[] den_mid;
delete[] velX_mid;
delete[] velY_mid;
delete[] flxX_mid;
delete[] flxY_mid;
delete[] vel_snd_arr;
delete[] vel_snd_arr_mid;
delete[] hden_mid;
//delete[] hvelX_mid;                    <= ??
//delete[] hvelY_mid;
delete[] hflxX_mid;
delete[] hflxY_mid;
delete[] hvel_snd_arr;
delete[] hvel_snd_arr_mid;

}

float DiracGraphene2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}

void DiracGraphene2D::ChooseGridPointers(const string &grid) {
	if(grid == "MidGrid"){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		ptr_snd = vel_snd_arr;
		ptr_den = Den;
		ptr_px = FlxX;
		ptr_py = FlxY;
		hptr_den = HDen;
		hptr_px = HFlxX;
		hptr_py = HFlxY;
		
	}if(grid == "MainGrid"){ // e vice-versa
		ptr_snd = vel_snd_arr_mid;
		ptr_den = den_mid;
		ptr_px = flxX_mid;
		ptr_py = flxY_mid;
		hptr_den = hden_mid;
		hptr_px = hflxX_mid;
		hptr_py = hflxY_mid;
		
	}
}


void DiracGraphene2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Den[c] = 1.0f + 0.005f * (noise - 0.5f);

		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        HDen[c] = 1.0f + 0.005f * (noise - 0.5f);
	}
}

void DiracGraphene2D::Richtmyer(){

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			GridPoint2D midpoint(ks, Nx, Ny, true);

			float den_avg   = 0.25f * (Den[midpoint.SW] + Den[midpoint.SE] + Den[midpoint.NW] + Den[midpoint.NE]);
			float flx_x_avg = 0.25f * (FlxX[midpoint.SW] + FlxX[midpoint.SE] + FlxX[midpoint.NW] + FlxX[midpoint.NE]);
			float flx_y_avg = 0.25f * (FlxY[midpoint.SW] + FlxY[midpoint.SE] + FlxY[midpoint.NW] + FlxY[midpoint.NE]);
			
            den_mid[ks] = den_avg
			              -0.5f*(dt/dx)*(DensityFluxX(midpoint, 'E') - DensityFluxX(midpoint, 'W'))
			              -0.5f*(dt/dy)*(DensityFluxY(midpoint, 'N') - DensityFluxY(midpoint, 'S'))
						  +0.5f*dt* DensitySource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
			flxX_mid[ks] = flx_x_avg
					-0.5f*(dt/dx)*(XMomentumFluxX(midpoint, 'E') - XMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(XMomentumFluxY(midpoint, 'N') - XMomentumFluxY(midpoint, 'S'))
					+0.5f*dt*XMomentumSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
			flxY_mid[ks] = flx_y_avg
					-0.5f*(dt/dx)*(YMomentumFluxX(midpoint, 'E') - YMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(YMomentumFluxY(midpoint, 'N') - YMomentumFluxY(midpoint, 'S'))
					+0.5f*dt*YMomentumSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);

			//____________________ HOLES _________________________
			
			float hden_avg   = 0.25f * (HDen[midpoint.SW] + HDen[midpoint.SE] + HDen[midpoint.NW] + HDen[midpoint.NE]);
			float hflx_x_avg = 0.25f * (HFlxX[midpoint.SW] + HFlxX[midpoint.SE] + HFlxX[midpoint.NW] + HFlxX[midpoint.NE]);
			float hflx_y_avg = 0.25f * (HFlxY[midpoint.SW] + HFlxY[midpoint.SE] + HFlxY[midpoint.NW] + HFlxY[midpoint.NE]);
			
            hden_mid[ks] = hden_avg
			              -0.5f*(dt/dx)*(HDensityFluxX(midpoint, 'E') - HDensityFluxX(midpoint, 'W'))
			              -0.5f*(dt/dy)*(HDensityFluxY(midpoint, 'N') - HDensityFluxY(midpoint, 'S'));
						  //+0.5f*dt* DensitySource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
			hflxX_mid[ks] = hflx_x_avg
					-0.5f*(dt/dx)*(HXMomentumFluxX(midpoint, 'E') - HXMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(HXMomentumFluxY(midpoint, 'N') - HXMomentumFluxY(midpoint, 'S'));
					//+0.5f*dt*XMomentumSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
			hflxY_mid[ks] = hflx_y_avg
					-0.5f*(dt/dx)*(HYMomentumFluxX(midpoint, 'E') - HYMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(HYMomentumFluxY(midpoint, 'N') - HYMomentumFluxY(midpoint, 'S'));
					//+0.5f*dt*YMomentumSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
		}
	
	ChooseGridPointers("MainGrid");
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			GridPoint2D mainpoint(kp, Nx, Ny, false);
			if( kp%Nx!=Nx-1 && kp%Nx!=0){
				float den_old = Den[kp];
				float flx_x_old = FlxX[kp];
				float flx_y_old = FlxY[kp];

                Den[kp] = den_old - (dt/dx)*(DensityFluxX(mainpoint, 'E') - DensityFluxX(mainpoint, 'W'))
						          - (dt/dy)*(DensityFluxY(mainpoint, 'N') - DensityFluxY(mainpoint, 'S'))
						          + dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxX[kp] = flx_x_old - (dt/dx)*(XMomentumFluxX(mainpoint, 'E') - XMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(XMomentumFluxY(mainpoint, 'N') - XMomentumFluxY(mainpoint, 'S'))
						             + dt*XMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxY[kp] = flx_y_old - (dt/dx)*(YMomentumFluxX(mainpoint, 'E') - YMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(YMomentumFluxY(mainpoint, 'N') - YMomentumFluxY(mainpoint, 'S'))
				                     + dt*YMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);

				//____________________ HOLES _________________________

				float hden_old = HDen[kp];
				float hflx_x_old = HFlxX[kp];
				float hflx_y_old = HFlxY[kp];

                HDen[kp] = hden_old - (dt/dx)*(HDensityFluxX(mainpoint, 'E') - HDensityFluxX(mainpoint, 'W'))
						          - (dt/dy)*(HDensityFluxY(mainpoint, 'N') - HDensityFluxY(mainpoint, 'S'));
						         // + dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				HFlxX[kp] = hflx_x_old - (dt/dx)*(HXMomentumFluxX(mainpoint, 'E') - HXMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(HXMomentumFluxY(mainpoint, 'N') - HXMomentumFluxY(mainpoint, 'S'));
						            // + dt*XMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				HFlxY[kp] = hflx_y_old - (dt/dx)*(HYMomentumFluxX(mainpoint, 'E') - HYMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(HYMomentumFluxY(mainpoint, 'N') - HYMomentumFluxY(mainpoint, 'S'));
				                    // + dt*YMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);


			}
		}
}
