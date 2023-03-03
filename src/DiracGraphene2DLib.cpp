/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/DiracGraphene2DLib.h"



DiracGraphene2D::DiracGraphene2D(SetUpParametersCNP &input_parameters) : Fluid2D(input_parameters) {

	hsize_t dimsf[2];
	dimsf[0] = static_cast<hsize_t>(Ny);
	dimsf[1] = static_cast<hsize_t>(Nx);
	DataspaceHDen = new DataSpace(RANK, dimsf );
	DataspaceHVelX = new DataSpace(RANK, dimsf );
	DataspaceHVelY = new DataSpace(RANK, dimsf );
	GrpHDen = nullptr;
	GrpHVelX = nullptr;
	GrpHVelY = nullptr;



	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency
	therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity

	vel_therm = input_parameters.ThermalVelocity ;
	A = input_parameters.Diffusive_sourceterm ;
	B = input_parameters.Creation_sourceterm ;

	char buffer [100];
	sprintf (buffer, "S=%.2fvF=%.2fvT=%.2fA=%.3fB=%.3fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, vel_therm, A, B, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
	file_infix = buffer;

	// main grid variables Nx*Ny
	HDen 		= new float[Nx * Ny]();
	HVelX 		= new float[Nx * Ny]();
	HVelY 		= new float[Nx * Ny]();

	HoleUmain = new StateVec2D[Nx * Ny]();
	HoleUmid = new StateVec2D[(Nx - 1) * (Ny - 1)]();
}

void DiracGraphene2D::SetSimulationTime(){
	this->SetTmax(5.0f+0.02f*vel_therm+20.0f/vel_therm);
}

void DiracGraphene2D::CflCondition(){ // Eventual redefinition

	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	float lambda;
	if(vel_therm<vel_snd){
		lambda=0.75f+sqrt(1.6f+3.0f*vel_snd*vel_snd);
	}else{
		lambda=0.75f+sqrt(1.6f+3.0f*vel_therm*vel_therm);
	}
	float dt1 = 0.5f * dx/abs(lambda);
	dt = dt1;

	
	//  CFL condition for (1,9) Weighted explicit method
	
	if(kin_vis>0.0){
		float dt2 = 0.8f*0.5f*dx*dx/kin_vis;
		dt = min(dt1, dt2);
	}
	
}



float DiracGraphene2D::EleDensityFluxX(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes) {
	float den=Uelec.n();
	float px=Uelec.px();
	return px / sqrt(den);
}

float DiracGraphene2D::EleDensityFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes) {
	float den=Uelec.n();
	float py=Uelec.py();
	return py / sqrt(den);
}

float DiracGraphene2D::EleXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes) {

	float sound=vel_snd;
	float den=Uelec.n();
	float hden=Uholes.n();
	float px=Uelec.px();
	float mass=DensityToMass(den);

	return px * px / mass + sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden);
}

float DiracGraphene2D::EleXMomentumFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes) {
	float den=Uelec.n();
	float px=Uelec.px();
	float py=Uelec.py();
	float mass=DensityToMass(den);
	return px * py / mass;
}


float DiracGraphene2D::EleYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes) {
	float sound=vel_snd;
	float den=Uelec.n();
	float hden=Uholes.n();
	float py=Uelec.py();
	float mass=DensityToMass(den);
	return py * py / mass + sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden);
}

float DiracGraphene2D::EleYMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes) {
	float den=Uelec.n();
	float px=Uelec.px();
	float py=Uelec.py();
	float mass=DensityToMass(den);

	return px * py / mass;
}


float DiracGraphene2D::EleDensitySource(StateVec2D Uelec , StateVec2D Uholes) {
	return A*(1.0f - Uelec.n()) + B*(Uelec.n()*Uholes.n() - 1.0f);
}

float DiracGraphene2D::EleXMomentumSource(StateVec2D Uelec , StateVec2D Uholes) {
	return 0.0f;
}

float DiracGraphene2D::EleYMomentumSource(StateVec2D Uelec , StateVec2D Uholes) {
	return 0.0f;
}


//________________________________ HOLE ____________________________

float DiracGraphene2D::HolDensityFluxX(StateVec2D Uelec , StateVec2D Uholes) {
	float hden=Uholes.n();
	float px=Uholes.px();
	return px / sqrt(hden);
}

float DiracGraphene2D::HolDensityFluxY(StateVec2D Uelec , StateVec2D Uholes) {
	float hden=Uholes.n();
	float py=Uholes.py();
	return py / sqrt(hden);
}

float DiracGraphene2D::HolXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes) {

	float sound=vel_snd;
	float den=Uelec.n();
	float hden=Uholes.n();
	float px=Uholes.px();
	float mass=DensityToMass(hden);

	return px * px / mass - sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden);
}

float DiracGraphene2D::HolXMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes) {
	float hden=Uholes.n();
	float px=Uholes.px();
	float py=Uholes.py();
	float mass=DensityToMass(hden);
	return px * py / mass;
}


float DiracGraphene2D::HolYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes) {
	float sound=vel_snd;
	float den=Uelec.n();
	float hden=Uholes.n();
	float py=Uholes.py();
	float mass=DensityToMass(hden);
	return py * py / mass - sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden);
}

float DiracGraphene2D::HolYMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes) {
	float hden=Uholes.n();
	float px=Uholes.px();
	float py=Uholes.py();
	float mass=DensityToMass(hden);
	return px * py / mass;
}

float DiracGraphene2D::HolDensitySource(StateVec2D Uelec , StateVec2D Uholes) {
	return A*(1.0f - Uholes.n()) + B*(Uelec.n()*Uholes.n() - 1.0f);
}

float DiracGraphene2D::HolXMomentumSource(StateVec2D Uelec , StateVec2D Uholes) {
	return 0.0f;
}

float DiracGraphene2D::HolYMomentumSource(StateVec2D Uelec , StateVec2D Uholes) {
	return 0.0f;
}


DiracGraphene2D::~DiracGraphene2D(){
delete[] Den;
delete[] VelX;
delete[] VelY;
delete[] CurX;
delete[] CurY;
delete[] HDen;
delete[] HVelX;
delete[] HVelY;
delete[] vel_snd_arr;
delete[] HoleUmain;
delete[] HoleUmid;
}

float DiracGraphene2D::DensityToMass(float density) {
	return sqrt(density*density*density);
}

void DiracGraphene2D::ChooseGridPointers(const string &grid) {
	if(grid == "MidGrid"){
		ptr_StateVec = Umain;
		ptr_StateVecHole = HoleUmain;
	}if(grid == "MainGrid"){
		ptr_StateVec = Umid;
		ptr_StateVecHole = HoleUmid;
	}
}


void DiracGraphene2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		//Den[c] = 1.0f + 0.05f * (noise - 0.5f);
		Umain[c].n()=1.0f + 0.05f * (noise - 0.5f);
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        //HDen[c] = 1.0f + 0.05f * (noise - 0.5f);
		HoleUmain[c].n()=1.0f + 0.05f * (noise - 0.5f);
	}
}

void DiracGraphene2D::InitialCondPointDen(){
	for (int c = 0; c < Nx*Ny; c++ ){
		//Den[c] = 1.0f;
        //HDen[c] = 1.0f;
		Umain[c].n()=1.0f;
		HoleUmain[c].n()=1.0f;
	}
	for (int c = 3*Ny/7; c < 4*Ny/7; c++){
		for (int g = 2*Nx/7; g < 3*Nx/7; g++){
			//Den[c*Nx+g] = 1.2f;
			//HDen[c*Nx+g] = 0.8f;
			Umain[c*Nx+g].n()=1.2f;
			HoleUmain[c*Nx+g].n()=0.8f;
		}
	}
	for (int c = 3*Ny/7; c < 4*Ny/7; c++){
		for (int g = 4*Nx/7; g < 5*Nx/7; g++){
			//Den[c*Nx+g] = 0.8f;
			//HDen[c*Nx+g] = 1.2f;
			Umain[c*Nx+g].n()=0.8f;
			HoleUmain[c*Nx+g].n()=1.2f;
		}
	}
}

void DiracGraphene2D::InitialCondUniform(){
	for (int c = 0; c < Nx*Ny; c++ ){
		//Den[c] = 1.0f;
        //HDen[c] = 1.0f;
		Umain[c].n()=1.0f;
		HoleUmain[c].n()=1.0f;
	}
}

void DiracGraphene2D::WriteFluidFile(float t){
	int j=Ny/2;
	int pos_end = Nx - 1 + j*Nx ;
	int pos_ini = j*Nx ;
	if(!isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].px()) || !isfinite(Umain[pos_end].px())){
			cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
			CloseHdf5File();
			exit(EXIT_FAILURE);
		}
	data_preview << t <<"\t"<< Umain[pos_ini] <<"\t"<<Umain[pos_end]<< "\t"<< HoleUmain[pos_ini] <<"\t"<< HoleUmain[pos_end]<< "\n";
}

/*
void DiracGraphene2D::Richtmyer(){

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,ptr_den,ptr_px,ptr_py,HDen,HFlxX,HFlxY,hptr_den,hptr_px,hptr_py,ptr_snd)
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			GridPoint2D midpoint(ks, Nx, Ny, true);

			//________ Electrons
			float den_avg   = 0.25f * (Den[midpoint.SW] + Den[midpoint.SE] + Den[midpoint.NW] + Den[midpoint.NE]);
			float flx_x_avg = 0.25f * (FlxX[midpoint.SW] + FlxX[midpoint.SE] + FlxX[midpoint.NW] + FlxX[midpoint.NE]);
			float flx_y_avg = 0.25f * (FlxY[midpoint.SW] + FlxY[midpoint.SE] + FlxY[midpoint.NW] + FlxY[midpoint.NE]);
			
			//________ Holes
			float hden_avg   = 0.25f * (HDen[midpoint.SW] + HDen[midpoint.SE] + HDen[midpoint.NW] + HDen[midpoint.NE]);
			float hflx_x_avg = 0.25f * (HFlxX[midpoint.SW] + HFlxX[midpoint.SE] + HFlxX[midpoint.NW] + HFlxX[midpoint.NE]);
			float hflx_y_avg = 0.25f * (HFlxY[midpoint.SW] + HFlxY[midpoint.SE] + HFlxY[midpoint.NW] + HFlxY[midpoint.NE]);

			//____________________ Electrons _________________________

            den_mid[ks] = den_avg
			              -0.5f*(dt/dx)*(EleDensityFluxX(midpoint, 'E', StateVec2D()) - EleDensityFluxX(midpoint, 'W', StateVec2D()))
			              -0.5f*(dt/dy)*(EleDensityFluxY(midpoint, 'N') - EleDensityFluxY(midpoint, 'S'));
			+0.5f * dt * EleDensitySource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);
			flxX_mid[ks] = flx_x_avg
					-0.5f*(dt/dx)*(EleXMomentumFluxX(midpoint, 'E') - EleXMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(EleXMomentumFluxY(midpoint, 'N') - EleXMomentumFluxY(midpoint, 'S'));
			+0.5f * dt * EleXMomentumSource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);
			flxY_mid[ks] = flx_y_avg
					-0.5f*(dt/dx)*(EleYMomentumFluxX(midpoint, 'E') - EleYMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(EleYMomentumFluxY(midpoint, 'N') - EleYMomentumFluxY(midpoint, 'S'));
			+0.5f * dt * EleYMomentumSource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);

			//____________________ HOLES _________________________
			
            hden_mid[ks] = hden_avg
			              -0.5f*(dt/dx)*(HolDensityFluxX(midpoint, 'E') - HolDensityFluxX(midpoint, 'W'))
			              -0.5f*(dt/dy)*(HolDensityFluxY(midpoint, 'N') - HolDensityFluxY(midpoint, 'S'));
			+0.5f * dt * HolDensitySource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);
			hflxX_mid[ks] = hflx_x_avg
					-0.5f*(dt/dx)*(HolXMomentumFluxX(midpoint, 'E') - HolXMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(HolXMomentumFluxY(midpoint, 'N') - HolXMomentumFluxY(midpoint, 'S'));
			+0.5f * dt * HolXMomentumSource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);
			hflxY_mid[ks] = hflx_y_avg
					-0.5f*(dt/dx)*(HolYMomentumFluxX(midpoint, 'E') - HolYMomentumFluxX(midpoint, 'W'))
					-0.5f*(dt/dy)*(HolYMomentumFluxY(midpoint, 'N') - HolYMomentumFluxY(midpoint, 'S'));
			+0.5f * dt * HolYMomentumSource(den_avg, flx_x_avg, flx_y_avg, hden_avg, hflx_x_avg, hflx_y_avg, 0.0f, 0.0f);
		}
	
	ChooseGridPointers("MainGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,ptr_den,ptr_px,ptr_py,HDen,HFlxX,HFlxY,hptr_den,hptr_px,hptr_py,ptr_snd)
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			GridPoint2D mainpoint(kp, Nx, Ny, false);
			if( kp%Nx!=Nx-1 && kp%Nx!=0){
				
				//________ Electrons
				float den_old = Den[kp];
				float flx_x_old = FlxX[kp];
				float flx_y_old = FlxY[kp];

				//________ Holes
				float hden_old = HDen[kp];
				float hflx_x_old = HFlxX[kp];
				float hflx_y_old = HFlxY[kp];

				//____________________ Electrons _________________________

                Den[kp] = den_old - (dt/dx)*(EleDensityFluxX(mainpoint, 'E', StateVec2D()) - EleDensityFluxX(mainpoint, 'W',
                                                                                                             StateVec2D()))
						          - (dt/dy)*(EleDensityFluxY(mainpoint, 'N') - EleDensityFluxY(mainpoint, 'S'));
				+ dt * EleDensitySource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);
				FlxX[kp] = flx_x_old - (dt/dx)*(EleXMomentumFluxX(mainpoint, 'E') - EleXMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(EleXMomentumFluxY(mainpoint, 'N') - EleXMomentumFluxY(mainpoint, 'S'));
				+ dt * EleXMomentumSource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);
				FlxY[kp] = flx_y_old - (dt/dx)*(EleYMomentumFluxX(mainpoint, 'E') - EleYMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(EleYMomentumFluxY(mainpoint, 'N') - EleYMomentumFluxY(mainpoint, 'S'));
				+ dt * EleYMomentumSource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);

				//____________________ HOLES _________________________

                HDen[kp] = hden_old - (dt/dx)*(HolDensityFluxX(mainpoint, 'E') - HolDensityFluxX(mainpoint, 'W'))
						          - (dt/dy)*(HolDensityFluxY(mainpoint, 'N') - HolDensityFluxY(mainpoint, 'S'));
				+ dt * HolDensitySource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);
				HFlxX[kp] = hflx_x_old - (dt/dx)*(HolXMomentumFluxX(mainpoint, 'E') - HolXMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(HolXMomentumFluxY(mainpoint, 'N') - HolXMomentumFluxY(mainpoint, 'S'));
				+ dt * HolXMomentumSource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);
				HFlxY[kp] = hflx_y_old - (dt/dx)*(HolYMomentumFluxX(mainpoint, 'E') - HolYMomentumFluxX(mainpoint, 'W'))
						        	 - (dt/dy)*(HolYMomentumFluxY(mainpoint, 'N') - HolYMomentumFluxY(mainpoint, 'S'));
				+ dt * HolYMomentumSource(den_old, flx_x_old, flx_y_old, hden_old, hflx_x_old, hflx_y_old, 0.0f, 0.0f);


			}
		}
}
*/

void DiracGraphene2D::RichtmyerStep1(){

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
		GridPoint2D midpoint(ks, Nx, Ny, true);

		//electrons
		StateVec2D eUavg{};
		eUavg = 0.25f * (Umain[midpoint.SW] + Umain[midpoint.SE] + Umain[midpoint.NW] + Umain[midpoint.NE]);
		StateVec2D eUNorth{};
		StateVec2D eUSouth{};
		StateVec2D eUEast{};
		StateVec2D eUWest{};
		eUNorth = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.NW]);
		eUSouth = 0.5f*(Umain[midpoint.SE]+Umain[midpoint.SW]);
		eUEast = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.SE]);
		eUWest = 0.5f*(Umain[midpoint.NW]+Umain[midpoint.SW]);
		//holes
		StateVec2D hUavg{};
		hUavg = 0.25f * (HoleUmain[midpoint.SW] + HoleUmain[midpoint.SE] + HoleUmain[midpoint.NW] + HoleUmain[midpoint.NE]);
		StateVec2D hUNorth{};
		StateVec2D hUSouth{};
		StateVec2D hUEast{};
		StateVec2D hUWest{};
		hUNorth = 0.5f*(HoleUmain[midpoint.NE]+HoleUmain[midpoint.NW]);
		hUSouth = 0.5f*(HoleUmain[midpoint.SE]+HoleUmain[midpoint.SW]);
		hUEast = 0.5f*(HoleUmain[midpoint.NE]+HoleUmain[midpoint.SE]);
		hUWest = 0.5f*(HoleUmain[midpoint.NW]+HoleUmain[midpoint.SW]);

		//____________________ Electrons _________________________

		Umid[ks].n() =  eUavg.n()
		                -0.5f*(dt/dx)*(EleDensityFluxX(eUEast,hUEast) - EleDensityFluxX(eUWest,hUWest))
		                -0.5f*(dt/dy)*(EleDensityFluxY(eUNorth,hUNorth) - EleDensityFluxY(eUSouth,hUSouth))+0.5f*dt* EleDensitySource(eUavg,hUavg);

		Umid[ks].px() = eUavg.px()
		                -0.5f*(dt/dx)*(EleXMomentumFluxX(eUEast,hUEast) - EleXMomentumFluxX(eUWest,hUWest))
		                -0.5f*(dt/dy)*(EleXMomentumFluxY(eUNorth,hUNorth) - EleXMomentumFluxY(eUSouth,hUSouth))+0.5f*dt*EleXMomentumSource(eUavg,hUavg);

		Umid[ks].py() = eUavg.py()
		                -0.5f*(dt/dx)*(EleYMomentumFluxX(eUEast,hUEast) - EleYMomentumFluxX(eUWest,hUWest))
		                -0.5f*(dt/dy)*(EleYMomentumFluxY(eUNorth,hUNorth) - EleYMomentumFluxY(eUSouth,hUSouth))+0.5f*dt*EleYMomentumSource(eUavg,hUavg);

		//________________________________ HOLE ____________________________

		HoleUmid[ks].n() =  hUavg.n()
		                    -0.5f*(dt/dx)*(HolDensityFluxX(eUEast,hUEast) - HolDensityFluxX(eUWest,hUWest))
		                    -0.5f*(dt/dy)*(HolDensityFluxY(eUNorth,hUNorth) - HolDensityFluxY(eUSouth,hUSouth))+0.5f*dt* HolDensitySource(eUavg,hUavg);

		HoleUmid[ks].px() = hUavg.px()
		                    -0.5f*(dt/dx)*(HolXMomentumFluxX(eUEast,hUEast) - HolXMomentumFluxX(eUWest,hUWest))
		                    -0.5f*(dt/dy)*(HolXMomentumFluxY(eUNorth,hUNorth) - HolXMomentumFluxY(eUSouth,hUSouth))+0.5f*dt*HolXMomentumSource(eUavg,hUavg);

		HoleUmid[ks].py() = hUavg.py()
		                    -0.5f*(dt/dx)*(HolYMomentumFluxX(eUEast,hUEast) - HolYMomentumFluxX(eUWest,hUWest))
		                    -0.5f*(dt/dy)*(HolYMomentumFluxY(eUNorth,hUNorth) - HolYMomentumFluxY(eUSouth,hUSouth))+0.5f*dt*HolYMomentumSource(eUavg,hUavg);
	}
}

void DiracGraphene2D::RichtmyerStep2(){

	ChooseGridPointers("MainGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		GridPoint2D mainpoint(kp, Nx, Ny, false);
		if( kp%Nx!=Nx-1 && kp%Nx!=0){

			//electrons
			StateVec2D eUold(Umain[kp]);
			StateVec2D eUNorth{};
			StateVec2D eUSouth{};
			StateVec2D eUEast{};
			StateVec2D eUWest{};
			eUNorth = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.NW]);
			eUSouth = 0.5f*(Umid[mainpoint.SE]+Umid[mainpoint.SW]);
			eUEast = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.SE]);
			eUWest = 0.5f*(Umid[mainpoint.NW]+Umid[mainpoint.SW]);

			//holes
			StateVec2D hUold(HoleUmain[kp]);
			StateVec2D hUNorth{};
			StateVec2D hUSouth{};
			StateVec2D hUEast{};
			StateVec2D hUWest{};
			hUNorth = 0.5f*(HoleUmid[mainpoint.NE]+HoleUmid[mainpoint.NW]);
			hUSouth = 0.5f*(HoleUmid[mainpoint.SE]+HoleUmid[mainpoint.SW]);
			hUEast = 0.5f*(HoleUmid[mainpoint.NE]+HoleUmid[mainpoint.SE]);
			hUWest = 0.5f*(HoleUmid[mainpoint.NW]+HoleUmid[mainpoint.SW]);


			//____________________ Electrons _________________________


			Umain[kp].n() = eUold.n()
			                - (dt/dx)*(EleDensityFluxX(eUEast,hUEast) - EleDensityFluxX(eUWest,hUWest))
			                - (dt/dy)*(EleDensityFluxY(eUNorth,hUNorth) - EleDensityFluxY(eUSouth,hUSouth))+ dt*EleDensitySource(eUold,hUold);
			Umain[kp].px() = eUold.px()
			                 - (dt/dx)*(EleXMomentumFluxX(eUEast,hUEast) - EleXMomentumFluxX(eUWest,hUWest))
			                 - (dt/dy)*(EleXMomentumFluxY(eUNorth,hUNorth) - EleXMomentumFluxY(eUSouth,hUSouth))+ dt*EleXMomentumSource(eUold,hUold);

			Umain[kp].py() = eUold.py()
			                 - (dt/dx)*(EleYMomentumFluxX(eUEast,hUEast) - EleYMomentumFluxX(eUWest,hUWest))
			                 - (dt/dy)*(EleYMomentumFluxY(eUNorth,hUNorth) - EleYMomentumFluxY(eUSouth,hUSouth))+ dt*EleYMomentumSource(eUold,hUold);


			//________________________________ HOLE ____________________________

			HoleUmain[kp].n() = hUold.n()
			                    - (dt/dx)*(HolDensityFluxX(eUEast,hUEast) - HolDensityFluxX(eUWest,hUWest))
			                    - (dt/dy)*(HolDensityFluxY(eUNorth,hUNorth) - HolDensityFluxY(eUSouth,hUSouth))+ dt*HolDensitySource(eUold,hUold);
			HoleUmain[kp].px() = hUold.px()
			                     - (dt/dx)*(HolXMomentumFluxX(eUEast,hUEast) - HolXMomentumFluxX(eUWest,hUWest))
			                     - (dt/dy)*(HolXMomentumFluxY(eUNorth,hUNorth) - HolXMomentumFluxY(eUSouth,hUSouth))+ dt*HolXMomentumSource(eUold,hUold);

			HoleUmain[kp].py() = hUold.py()
			                     - (dt/dx)*(HolYMomentumFluxX(eUEast,hUEast) - HolYMomentumFluxX(eUWest,hUWest))
			                     - (dt/dy)*(HolYMomentumFluxY(eUNorth,hUNorth) - HolYMomentumFluxY(eUSouth,hUSouth))+ dt*HolYMomentumSource(eUold,hUold);
		}
	}
}

void DiracGraphene2D::ForwardTimeOperator() {
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old,hflx_x_old, hflx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			flx_x_old = Umain[kp].px();
			flx_y_old = Umain[kp].py();
			Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
			Umain[kp].py() = flx_y_old + Umain[kp].d2vy();

			hflx_x_old = HoleUmain[kp].px();
			hflx_y_old = HoleUmain[kp].py();
			HoleUmain[kp].px() = hflx_x_old + HoleUmain[kp].d2vx();
			HoleUmain[kp].py() = hflx_y_old + HoleUmain[kp].d2vy();

			tmp_old = Umain[kp].tmp();
			Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();


		}
	}
}
void DiracGraphene2D::ForwardTimeOperator(char field) {
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old,hflx_x_old, hflx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			switch(field) {
				case 'T': 	tmp_old = Umain[kp].tmp();
					Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();
					break;
				case 'V':
					flx_x_old = Umain[kp].px();
					flx_y_old = Umain[kp].py();
					Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
					Umain[kp].py() = flx_y_old + Umain[kp].d2vy();
					hflx_x_old = HoleUmain[kp].px();
					hflx_y_old = HoleUmain[kp].py();
					HoleUmain[kp].px() = hflx_x_old + HoleUmain[kp].d2vx();
					HoleUmain[kp].py() = hflx_y_old + HoleUmain[kp].d2vy();
					break;
				default: ;
			}
		}
	}
}



void DiracGraphene2D::VelocityLaplacianWeighted19() {


	for (int kp = 0; kp < Nx * Ny ; kp++) {
		float den =Umain[kp].n();
		float mass = DensityToMass(den);
		VelX[kp] = Umain[kp].px()/mass;
		VelY[kp] = Umain[kp].py()/mass;
		float hden =HoleUmain[kp].n();
		float hmass = DensityToMass(hden);
		HVelX[kp] = HoleUmain[kp].px()/hmass;
		HVelY[kp] = HoleUmain[kp].py()/hmass;
	}

#pragma omp parallel for default(none) shared(Umain,VelX,VelY)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		GridPoint2D point(kp,Nx,Ny,false);
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			Umain[kp].d2vx() = Laplacian19( point, VelX, kin_vis);
			Umain[kp].d2vy() = Laplacian19( point, VelY, kin_vis);
			HoleUmain[kp].d2vx() = Laplacian19( point, HVelX, kin_vis);
			HoleUmain[kp].d2vy() = Laplacian19( point, HVelY, kin_vis);
		}
	}
}

void DiracGraphene2D::CopyFields() {
	float mass, hmass;
	for (int i = 0; i < Nx*Ny; ++i) {
		Den[i]=Umain[i].n();
		mass= DensityToMass(Den[i]);
		VelX[i]=Umain[i].px()/mass;
		VelY[i]=Umain[i].py()/mass;
		Tmp[i] =Umain[i].tmp();
		HDen[i]=HoleUmain[i].n();
		hmass= DensityToMass(HDen[i]);
		HVelX[i]=HoleUmain[i].px()/hmass;
		HVelY[i]=HoleUmain[i].py()/hmass;
	}
}

void DiracGraphene2D::SaveSnapShot() {

	Fluid2D::SaveSnapShot();

	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;

	this->CopyFields();

	string str_time = to_string(TimeStepCounter / snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

	float currenttime=static_cast<float>(TimeStepCounter) * dt;

	DataSet dataset_hden = GrpHDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceHDen);
	Attribute atr_step_hden = dataset_hden.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_hden = dataset_hden.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_hden.write(HDen, HDF5FLOAT);
	dataset_hden.close();
	atr_step_hden.write(HDF5INT, &TimeStepCounter);
	atr_time_hden.write(HDF5FLOAT , &currenttime);
	atr_step_hden.close();
	atr_time_hden.close();


	DataSet dataset_hvel_x = GrpHVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceHVelX);
	Attribute atr_step_hvel_x = dataset_hvel_x.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_hvel_x = dataset_hvel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_hvel_x.write(VelX, HDF5FLOAT);
	dataset_hvel_x.close();
	atr_step_hvel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_hvel_x.write(HDF5FLOAT , &currenttime);
	atr_step_hvel_x.close();
	atr_time_hvel_x.close();

	DataSet dataset_hvel_y = GrpHVelY->createDataSet(name_dataset, HDF5FLOAT, *DataspaceHVelY);
	Attribute atr_step_hvel_y = dataset_hvel_y.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_hvel_y = dataset_hvel_y.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_hvel_y.write(VelY, HDF5FLOAT);
	dataset_hvel_y.close();
	atr_step_hvel_y.write(HDF5INT, &TimeStepCounter);
	atr_time_hvel_y.write(HDF5FLOAT , &currenttime);
	atr_step_hvel_y.close();
	atr_time_hvel_y.close();

}

void DiracGraphene2D::CreateHdf5File() {
	TethysBase::CreateHdf5File();
	GrpHDen = new Group(Hdf5File->createGroup("/Data/HoleDensity" ));
	GrpHVelX = new Group(Hdf5File->createGroup("/Data/HoleVelocityX" ));
	GrpHVelY = new Group(Hdf5File->createGroup("/Data/HoleVelocityY"));
}

