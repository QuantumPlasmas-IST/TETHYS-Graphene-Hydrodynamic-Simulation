/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/Fluid2DLib.h"
#include "includes/SetUpParametersLib.h"


using namespace H5;
using namespace std;

Fluid2D::Fluid2D(const SetUpParameters &input_parameters) : TethysBase{input_parameters.SizeX, input_parameters.SizeY, 2}{
	Nx = input_parameters.SizeX;
	Ny = input_parameters.SizeY;
	lengX=input_parameters.Length;
	lengY=input_parameters.Width;
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	vel_snd = input_parameters.SoundVelocity;//sound_velocity;
	kin_vis = input_parameters.ShearViscosity;//shear_viscosity;
	odd_vis = input_parameters.OddViscosity;//odd_viscosity;

	if(input_parameters.SimulationTime > 0.0f){
		Tmax = input_parameters.SimulationTime;
	}

	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fodd=%.3fl=%.2fwc=%.2f", vel_snd, vel_fer, kin_vis,odd_vis, col_freq,cyc_freq);
	file_infix = buffer;
	// main grid variables Nx*Ny
	Tmp 		= new float[Nx * Ny]();
	Den 		= new float[Nx * Ny]();
	VelX 		= new float[Nx * Ny]();
	VelY 		= new float[Nx * Ny]();
/*
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx * Ny]();
	velX_dx	= new float[Nx * Ny]();
	velX_dy	= new float[Nx * Ny]();
	velY_dx	= new float[Nx * Ny]();
	velY_dy	= new float[Nx * Ny]();
*/
	den_dx = new float[Nx * Ny]();
	den_dy = new float[Nx * Ny]();

	lap_den = new float[Nx*Ny]();
	lap_den_mid =new float[(Nx-1)*(Ny-1)]();

	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve
    lap_tmp = new float[Nx*Ny]();

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
/*	tmp_mid		= new float[(Nx-1)*(Ny-1)]();
	den_mid		= new float[(Nx-1)*(Ny-1)]();
	velX_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_mid	= new float[(Nx-1)*(Ny-1)]();
	velX_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velX_dy_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dy_mid	= new float[(Nx-1)*(Ny-1)]();
*/
	den_dx_mid = new float[(Nx-1)*(Ny-1)]();
	den_dy_mid = new float[(Nx-1)*(Ny-1)]();

//	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
//	flxY_mid	= new float[(Nx-1)*(Ny-1)]();
//	vel_snd_arr_mid	= new float[(Nx-1)*(Ny-1)]();


	Umain = new StateVec2D[Nx*Ny]();
	Umid = new StateVec2D[(Nx-1)*(Ny-1)]();

}

Fluid2D::~Fluid2D() = default;

void Fluid2D::SetSound(){
	for(int kp=0; kp<=Nx*Ny-1; kp++) {
//		vel_snd_arr[kp] = vel_snd;
		Umain[kp].S() = vel_snd;
	}
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) {
//		vel_snd_arr_mid[ks]= vel_snd;
		Umid[ks].S()= vel_snd;
	}
}

void Fluid2D::SetSound(std::function<float(float,float)> func){
	for(int kp=0; kp<=Nx*Ny-1; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		//vel_snd_arr[kp]= func(i*dx,j*dy);
		Umain[kp].S()= func(i*dx,j*dy);
	}
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) { //correr todos os pontos da grelha secundaria
		div_t divresult;
		divresult = div(ks, Nx - 1);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		//vel_snd_arr_mid[ks]=  func((i+0.5f)*dx,(j+0.5f)*dy);
		Umid[ks].S()=  func((i+0.5f)*dx,(j+0.5f)*dy);
	}
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Initial condition setting
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void Fluid2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Umain[c].n() = 1.0f + 0.005f * (noise - 0.5f);
	//	Den[c] = 1.0f + 0.005f * (noise - 0.5f);
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Umain[c].tmp()=  .2f + 0.005f * (noise - 0.5f);
	//	Tmp[c]=  .2f + 0.005f * (noise - 0.5f);
	}
}

void Fluid2D::InitialCondGeneral(function<float(float, float)> fden, function<float(float, float)> fvx, function<float(
		float, float)> fvy) {
	float x,y;
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
			x=i*dx;
			y=j*dy;
			Umain[i + j * Nx].n() = fden(x,y);
			Umain[i + j * Nx].px() = fvx(x,y);
			Umain[i + j * Nx].py() = fvy(x,y);
		}
	}
}


void Fluid2D::InitialCondWave() {
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
			//Den[i + j * Nx] = 1.0f+0.3f*sin(2.0f*MAT_PI*i*dx/lengX) ;
			Umain[i + j * Nx].n() = 1.0f+0.05f*sin(5.0f*2.0f*MAT_PI*i*dx/lengX) ;
			//Den[i + j * Nx] = 1.0f+0.3f/cosh((i*dx-0.5f)*20.0f);
			Umain[i + j * Nx].px() = 0.0f;
		}
	}
}

void Fluid2D::InitialCondTest(){
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
			float densi;
			if(i>=80&&i<=120&&j>=80&&j<=120){
			densi=0.2f;
			}
			else{
			densi=0.0f;
			}
			Umain[i + j * Nx].n() = 1.0f + densi;
			Umain[i + j * Nx].px() = 0.1f;
		//	Den[i + j * Nx] = 1.0f + densi;
		//	FlxX[i + j * Nx] = 0.1f;
		}
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
void Fluid2D::MassFluxToVelocity(const string grid) {
	float den;
	if(grid=="MainGrid"){
		for(int c=0; c <= Nx * Ny - 1; c++){
			den = Den[c];
			VelX[c]= FlxX[c] / DensityToMass(den);
			VelY[c]= FlxY[c] / DensityToMass(den);
		}
	}
	if(grid=="MidGrid"){
		for(int c=0; c <= (Nx-1) * (Ny-1) - 1; c++){
			den = den_mid[c];
			velX_mid[c]= flxX_mid[c] / DensityToMass(den);
			velY_mid[c]= flxY_mid[c] / DensityToMass(den);
		}
	}
}

void Fluid2D::VelocityToCurrent() {
	for(int c=0; c <= Nx * Ny - 1; c++){
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}
*/


void Fluid2D::Richtmyer(){
	/*
	if(odd_vis) {
		//this->VelocityGradient();
		this->MassFluxToVelocity("MainGrid");
		GradientField(VelX,velX_dx,velX_dy,dx,dy,Nx,Ny); // TODO é capaz de ser mais intelignente separar estes dois calculos de gradiente por metade dos cores ja que sao calculos independente e mais ou menos equivalentes
		GradientField(VelY,velY_dx,velY_dy,dx,dy,Nx,Ny);
	}
	if(therm_diff) {
		//this->DensityGradient();
		GradientField(Den,den_dx,den_dy,dx,dy,Nx,Ny);
	}
	//this->DensityLaplacian();
	LaplacianField(Den,lap_den,dx,Nx,Ny);
	*/

	RichtmyerStep1();

	/*
	if(odd_vis) {
		//this->VelocityGradientMid();
		this->MassFluxToVelocity("MidGrid");
		GradientField(velX_mid,velX_dx_mid,velX_dy_mid,dx,dy,Nx-1,Ny-1); // TODO é capaz de ser mais intelignente separar estes dois calculos de gradiente por metade dos cores ja que sao calculos independente e mais ou menos equivalentes
		GradientField(velY_mid,velY_dx_mid,velY_dy_mid,dx,dy,Nx-1,Ny-1);

	}
	if(therm_diff) {
//		this->DensityGradient();
	GradientField(den_mid,den_dx_mid,den_dy_mid,dx,dy,Nx-1,Ny-1);
	}
	LaplacianField(den_mid,lap_den_mid,dx,Nx-1,Ny-1);
	*/

	RichtmyerStep2();


}



void Fluid2D::RichtmyerStep1(){

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
		GridPoint2D midpoint(ks, Nx, Ny, true);

		/*
		float den_avg   = 0.25f * (Den[midpoint.SW] + Den[midpoint.SE] + Den[midpoint.NW] + Den[midpoint.NE]);
		float flx_x_avg = 0.25f * (FlxX[midpoint.SW] + FlxX[midpoint.SE] + FlxX[midpoint.NW] + FlxX[midpoint.NE]);
		float flx_y_avg = 0.25f * (FlxY[midpoint.SW] + FlxY[midpoint.SE] + FlxY[midpoint.NW] + FlxY[midpoint.NE]);
		float tmp_avg   = 0.25f * (Tmp[midpoint.SW] + Tmp[midpoint.SE] + Tmp[midpoint.NW] + Tmp[midpoint.NE]);
		*/



		/*
		float n_dx=0.0f,n_dy=0.0f;
		if(therm_diff) {
			n_dx =0.25f * (den_dx[midpoint.SW] + den_dx[midpoint.SE] + den_dx[midpoint.NW] + den_dx[midpoint.NE]);
			n_dy =0.25f * (den_dy[midpoint.SW] + den_dy[midpoint.SE] + den_dy[midpoint.NW] + den_dy[midpoint.NE]) ;
		}
		*/


		StateVec2D Uavg{};
		Uavg = 0.25f * (Umain[midpoint.SW] + Umain[midpoint.SE] + Umain[midpoint.NW] + Umain[midpoint.NE]);

		StateVec2D UNorth{};
		StateVec2D USouth{};
		StateVec2D UEast{};
		StateVec2D UWest{};
		UNorth = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.NW]);
		USouth = 0.5f*(Umain[midpoint.SE]+Umain[midpoint.SW]);
		UEast = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.SE]);
		UWest = 0.5f*(Umain[midpoint.NW]+Umain[midpoint.SW]);

/*		UEast = SideAverage(ptr_StateVec,midpoint,'E');
		UWest = SideAverage(ptr_StateVec,midpoint,'W');
		UNorth = SideAverage(ptr_StateVec,midpoint,'N');
		USouth = SideAverage(ptr_StateVec,midpoint,'S');
*/
		Umid[ks].n() =  Uavg.n()
		                -0.5f*(dt/dx)*(DensityFluxX(UEast) - DensityFluxX(UWest))
		                -0.5f*(dt/dy)*(DensityFluxY(UNorth) - DensityFluxY(USouth));
		                //+0.5f*dt* DensitySource(Uavg);

		Umid[ks].px() = Uavg.px()
		                -0.5f*(dt/dx)*(XMomentumFluxX(UEast) - XMomentumFluxX(UWest))
		                -0.5f*(dt/dy)*(XMomentumFluxY(UNorth) - XMomentumFluxY(USouth));
		                //+0.5f*dt*XMomentumSource(Uavg);

		Umid[ks].py() = Uavg.py()
		                -0.5f*(dt/dx)*(YMomentumFluxX(UEast) - YMomentumFluxX(UWest))
		                -0.5f*(dt/dy)*(YMomentumFluxY(UNorth) - YMomentumFluxY(USouth));
		                //+0.5f*dt*YMomentumSource(Uavg);

/*
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
*/



		/*
		if(therm_diff){
			tmp_mid[ks] = tmp_avg
			              -0.5f * (dt / dx) * (TemperatureFluxX(midpoint, 'E') - TemperatureFluxX(midpoint, 'W'))
			              -0.5f * (dt / dy) * (TemperatureFluxY(midpoint, 'N') - TemperatureFluxY(midpoint, 'S'))
			              +0.5f * dt * TemperatureSource(den_avg, flx_x_avg, flx_y_avg, n_dx, n_dy, 0.0f, 0.0f);
		}*/
	}



}


void Fluid2D::RichtmyerStep2(){

	ChooseGridPointers("MainGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		GridPoint2D mainpoint(kp, Nx, Ny, false);
		if( kp%Nx!=Nx-1 && kp%Nx!=0){


		/*	float den_old = Den[kp];
			float flx_x_old = FlxX[kp];
			float flx_y_old = FlxY[kp];
			float tmp_old = Tmp[kp];
		 */

			StateVec2D Uold(Umain[kp]);
			StateVec2D UNorth{};
			StateVec2D USouth{};
			StateVec2D UEast{};
			StateVec2D UWest{};

			UNorth = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.NW]);
			USouth = 0.5f*(Umid[mainpoint.SE]+Umid[mainpoint.SW]);
			UEast = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.SE]);
			UWest = 0.5f*(Umid[mainpoint.NW]+Umid[mainpoint.SW]);


	//		UEast = SideAverage(ptr_StateVec,mainpoint,'E');
	//		UWest = SideAverage(ptr_StateVec,mainpoint,'W');
	//		UNorth = SideAverage(ptr_StateVec,mainpoint,'N');
	//		USouth = SideAverage(ptr_StateVec,mainpoint,'S');

			Umain[kp].n() = Uold.n()
			                - (dt/dx)*(DensityFluxX(UEast) - DensityFluxX(UWest))
			                - (dt/dy)*(DensityFluxY(UNorth) - DensityFluxY(USouth));
			                //+ dt*DensitySource(Uold);
			Umain[kp].px() = Uold.px()
			                 - (dt/dx)*(XMomentumFluxX(UEast) - XMomentumFluxX(UWest))
			                 - (dt/dy)*(XMomentumFluxY(UNorth) - XMomentumFluxY(USouth));
			                 //+ dt*XMomentumSource(Uold);

			Umain[kp].py() = Uold.py()
			                 - (dt/dx)*(YMomentumFluxX(UEast) - YMomentumFluxX(UWest))
			                 - (dt/dy)*(YMomentumFluxY(UNorth) - YMomentumFluxY(USouth));
			                 //+ dt*YMomentumSource(Uold);







/*
			Den[kp] = den_old - (dt/dx)*(DensityFluxX(mainpoint, 'E') - DensityFluxX(mainpoint, 'W'))
			          - (dt/dy)*(DensityFluxY(mainpoint, 'N') - DensityFluxY(mainpoint, 'S'))
			          + dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
			FlxX[kp] = flx_x_old - (dt/dx)*(XMomentumFluxX(mainpoint, 'E') - XMomentumFluxX(mainpoint, 'W'))
			           - (dt/dy)*(XMomentumFluxY(mainpoint, 'N') - XMomentumFluxY(mainpoint, 'S'))
			           + dt*XMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
			FlxY[kp] = flx_y_old - (dt/dx)*(YMomentumFluxX(mainpoint, 'E') - YMomentumFluxX(mainpoint, 'W'))
			           - (dt/dy)*(YMomentumFluxY(mainpoint, 'N') - YMomentumFluxY(mainpoint, 'S'))
			           + dt*YMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);

*/


			/*
			if(therm_diff) {
				Tmp[kp] = tmp_old - (dt / dx) * (TemperatureFluxX(mainpoint, 'E') - TemperatureFluxX(mainpoint, 'W'))
				          - (dt / dy) * (TemperatureFluxY(mainpoint, 'N') - TemperatureFluxY(mainpoint, 'S'))
				          + dt * TemperatureSource(den_old, flx_x_old, flx_y_old, den_dx[kp], den_dy[kp], 0.0f, 0.0f);
			}
			*/
		}
	}

}



void Fluid2D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0f;
}

void Fluid2D::CreateFluidFile(){
	std::string previewfile = "preview_2D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid2D::WriteFluidFile(float t){
	int j=Ny/2;
	int pos_end = Nx - 1 + j*Nx ;
	int pos_ini = j*Nx ;
		//if(!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(FlxX[pos_end]) || !isfinite(FlxX[pos_ini])){
		if(!isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].px()) || !isfinite(Umain[pos_end].px())){
			cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
			CloseHdf5File();
			exit(EXIT_FAILURE);
		}

/*	data_preview << t << "\t"
	<< Den[pos_end]  << "\t"
	<< FlxX[pos_end] << "\t"
	<< Den[pos_ini]  << "\t"
	<< FlxX[pos_ini] << "\t\t";*/
	data_preview << t <<"\t"<< Umain[pos_ini] <<"\t"<<Umain[pos_end]<< "\n";
}

void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}

void Fluid2D::VelocityLaplacianFtcs() {
	this->MassFluxToVelocity("MainGrid");

//calculate laplacians
#pragma omp parallel for  default(none) shared(lap_flxX,lap_flxY,VelX,VelY,Nx,Ny)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		int north, south, east, west;
		div_t divresult;
		divresult = div(kp, Nx);
		int j;
		j = divresult.quot;
		int i;
		i = divresult.rem;
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			north = i + (j + 1) * Nx;
			south = i + (j - 1) * Nx;
			east = i + 1 + j * Nx;
			west = i - 1 + j * Nx;
			lap_flxX[kp] =
					kin_vis*dt*(-4.0f * VelX[kp]  + VelX[north]  + VelX[south]  + VelX[east]  +
							VelX[west] ) / (dx * dx);
			lap_flxY[kp] =
					kin_vis*dt*(-4.0f * VelY[kp]  + VelY[north]  + VelY[south]  + VelY[east]  +
							VelY[west] ) / (dx * dx);
		}
	}
}

void Fluid2D::VelocityLaplacianWeighted19() {
	this->MassFluxToVelocity("MainGrid");

#pragma omp parallel for default(none) shared(lap_flxX,lap_flxY,VelX,VelY)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		GridPoint2D point(kp,Nx,Ny,false);
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			lap_flxX[kp] = Laplacian19( point, VelX, kin_vis);
			lap_flxY[kp] = Laplacian19( point, VelY, kin_vis);
		}
	}
}

void Fluid2D::TemperatureLaplacianWeighted19() {
#pragma omp parallel for default(none) shared(lap_tmp,Tmp)
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
	    GridPoint2D point(kp,Nx,Ny,false);
        if (kp % Nx != Nx - 1 && kp % Nx != 0){
	        lap_tmp[kp] = Laplacian19( point, Tmp, therm_diff );
        }
    }
}


void Fluid2D:: ParabolicOperatorFtcs() {
	VelocityLaplacianFtcs();
	ForwardTimeOperator();
}

void Fluid2D::ParabolicOperatorWeightedExplicit19() {
	VelocityLaplacianWeighted19();
	TemperatureLaplacianWeighted19();
	ForwardTimeOperator();
}
void Fluid2D::ParabolicOperatorWeightedExplicit19(char field) {
	switch(field) {
		case 'V': VelocityLaplacianWeighted19();
			ForwardTimeOperator(field);
			break;
		case 'T': TemperatureLaplacianWeighted19();
			ForwardTimeOperator(field);
			break;
		default: ;
	}
}

void Fluid2D::SaveSound() {

	for (int i = 0; i < Nx*Ny; ++i) {
		vel_snd_arr[i]=Umain[i].S();
	}

	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocicity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(vel_snd_arr, HDF5FLOAT);
	dataset_vel_snd.close();
}

void Fluid2D::ReadSnapShot(const H5std_string &snap_name) {
	DataSet dataset_den = GrpDen->openDataSet( snap_name );
	DataSet dataset_vel_x  = GrpVelX->openDataSet( snap_name );
	DataSet dataset_vel_y = GrpVelY->openDataSet( snap_name );
	DataSpace dataspace = dataset_den.getSpace(); // Get dataspace of the dataset.
	if(dataset_den.attrExists ("time")){
		Attribute attr_time = dataset_den.openAttribute("time");
		attr_time.read(attr_time.getDataType(), &TimeStamp);
	}else{
		TimeStamp+=1.0f;
	}
	if(dataset_den.attrExists ("time step")){
		Attribute attr_counter = dataset_den.openAttribute("time step");
		attr_counter.read(attr_counter.getDataType(), &TimeStepCounter);
	}else{
		TimeStepCounter++;
	}
	dataset_den.read( Den, PredType::NATIVE_FLOAT,   dataspace );
	dataset_vel_x.read( VelX, PredType::NATIVE_FLOAT,   dataspace );
	dataset_vel_y.read( VelY, PredType::NATIVE_FLOAT,   dataspace );
	dataset_den.close();
	dataset_vel_x.close();
	dataset_vel_y.close();
}

void Fluid2D::SaveSnapShot() {
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;

	//this->MassFluxToVelocity("MainGrid");

	this->CopyFields();

	string str_time = to_string(TimeStepCounter / snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

    float currenttime=static_cast<float>(TimeStepCounter) * dt;

    DataSet dataset_den = GrpDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_den.write(Den, HDF5FLOAT);
    dataset_den.close();
	atr_step_den.write(HDF5INT, &TimeStepCounter);
	atr_time_den.write(HDF5FLOAT , &currenttime);
	atr_step_den.close();
	atr_time_den.close();



	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelX);
	Attribute atr_step_vel_x = dataset_vel_x.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_x = dataset_vel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_vel_x.write(VelX, HDF5FLOAT);
	dataset_vel_x.close();
	atr_step_vel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_x.write(HDF5FLOAT , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();

	DataSet dataset_vel_y = GrpVelY->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelY);
	Attribute atr_step_vel_y = dataset_vel_y.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_y = dataset_vel_y.createAttribute("time", HDF5FLOAT, atr_dataspace);
	dataset_vel_y.write(VelY, HDF5FLOAT);
	dataset_vel_y.close();
	atr_step_vel_y.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_y.write(HDF5FLOAT , &currenttime);
	atr_step_vel_y.close();
	atr_time_vel_y.close();

    DataSet dataset_tmp = GrpTmp->createDataSet(name_dataset, HDF5FLOAT, *DataspaceTmp);
    Attribute atr_step_tmp = dataset_tmp.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_tmp = dataset_tmp.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_tmp.write(Tmp, HDF5FLOAT);
    dataset_tmp.close();
    atr_step_tmp.write(HDF5INT, &TimeStepCounter);
    atr_time_tmp.write(HDF5FLOAT , &currenttime);
    atr_step_tmp.close();
    atr_time_tmp.close();
}

int Fluid2D::GetSnapshotStep() const { return snapshot_step;}
int Fluid2D::GetSnapshotFreq() const {return snapshot_per_period;}

bool Fluid2D::Snapshot() const {
	bool state;
	if(TimeStepCounter % snapshot_step == 0){
		state = true;
	}else{
		state = false;
	}
	return state;
}

float Fluid2D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

float Fluid2D::XMomentumSource(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return 0.0f;
}
float Fluid2D::YMomentumSource(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s) {
	return 0.0f;
}

void Fluid2D::ForwardTimeOperator() {
//#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			//flx_x_old = FlxX[kp];
			//flx_y_old = FlxY[kp];

			flx_x_old = Umain[kp].px();
			flx_y_old = Umain[kp].py();
            tmp_old = Umain[kp].tmp();

			Umain[kp].px() = flx_x_old + lap_flxX[kp];
			Umain[kp].py() = flx_y_old + lap_flxY[kp];
			Umain[kp].tmp() = tmp_old + lap_tmp[kp];
        }
	}
}
void Fluid2D::ForwardTimeOperator(char field) {  //TODO meter o switch
//#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			//flx_x_old = FlxX[kp];
			//flx_y_old = FlxY[kp];

			flx_x_old = Umain[kp].px();
			flx_y_old = Umain[kp].py();
			tmp_old = Umain[kp].tmp();

			Umain[kp].px() = flx_x_old + lap_flxX[kp];
			Umain[kp].py() = flx_y_old + lap_flxY[kp];
			Umain[kp].tmp() = tmp_old + lap_tmp[kp];
		}
	}
}

float Fluid2D::DensityToMass(float density) {
	return density;
}
/*
float Fluid2D::DensityFluxX(GridPoint2D p, char side) {
	float px;
	px = SideAverage(ptr_px, p, side);
	return px;
}*/
float Fluid2D::DensityFluxX(StateVec2D U) {
	return U.px();
}

/*
float Fluid2D::DensityFluxY(GridPoint2D p, char side) {
	float py;
	py = SideAverage(ptr_py, p, side);
	return py;
}*/

float Fluid2D::DensityFluxY(StateVec2D U) {
	return U.py();
}

/*
float Fluid2D::XMomentumFluxX(GridPoint2D p, char side) {
	float den;
	float px;
	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	return px * px / den + den;
}*/
float Fluid2D::XMomentumFluxX(StateVec2D U) {
	return U.px()*U.px()/U.n() + U.n();
}
/*
float Fluid2D::XMomentumFluxY(GridPoint2D p, char side) {
	float den;
	float px;
	float py;
	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	return px * py / den;
}*/

float Fluid2D::XMomentumFluxY(StateVec2D U) {
	return U.px()*U.py()/U.n();
}


/*
float Fluid2D::YMomentumFluxY(GridPoint2D p, char side) {
	float den;
	float py;
	den = SideAverage(ptr_den, p, side);
	py = SideAverage(ptr_py, p, side);
	return py * py / den + den;
}*/
float Fluid2D::YMomentumFluxY(StateVec2D U) {
	return U.py()*U.py()/U.n() + U.n();
}

/*
float Fluid2D::YMomentumFluxX(GridPoint2D p, char side) {
	float den;
	float px;
	float py;
	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	py = SideAverage(ptr_py, p, side);
	return px * py / den;
}*/

float Fluid2D::YMomentumFluxX(StateVec2D U) {
	return U.py()*U.px()/U.n();
}

/*
float Fluid2D::TemperatureFluxX(GridPoint2D p, char side) {
	float tmp;
    float px;
	float den;

	den = SideAverage(ptr_den, p, side);
	px = SideAverage(ptr_px, p, side);
	tmp = SideAverage(ptr_tmp, p, side);
	return PHYS_FERMI_CNVC*px  + tmp*px / DensityToMass(den);
}


float Fluid2D::TemperatureFluxY(GridPoint2D p, char side) {
    float py;
	float tmp;
	float den;

	den = SideAverage(ptr_den, p, side);
	py = SideAverage(ptr_py, p, side);
	tmp = SideAverage(ptr_tmp, p, side);
	return PHYS_FERMI_CNVC*py  + tmp*py / DensityToMass(den) ;
}

*/

float Fluid2D::Laplacian19(GridPoint2D p, float *input_ptr, float constant) {
	float sx=constant*dt/(dx*dx);
	float sy=constant*dt/(dy*dy);
	float * data_ptr = input_ptr;
	float lap;
	lap = (4.0f*sx*sy-2.0f*sx-2.0f*sy)*data_ptr[p.C] +
	               sx*sy*( data_ptr[p.NE2] + data_ptr[p.SE2] + data_ptr[p.NW2] + data_ptr[p.SW2])
	               + sy*(1.0f-2.0f*sx)*(data_ptr[p.N] + data_ptr[p.S])
	               + sx*(1.0f-2.0f*sy)*(data_ptr[p.W] + data_ptr[p.E]);
return lap;
}

float Fluid2D::TemperatureSource(float n, float flx_x, float flx_y, float den_grad_x, float den_grad_y, float mass, float s) {
	return 0;
}

void Fluid2D::ChooseGridPointers(const string &grid) {
	if(grid == "MidGrid"){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
/*		ptr_snd = vel_snd_arr;
		ptr_den = Den;
		ptr_px = FlxX;
		ptr_py = FlxY;
		ptr_dendx = den_dx;
		ptr_dendy = den_dy;
		ptr_velXdx = velX_dx;
		ptr_velXdy = velX_dy;
		ptr_velYdx = velY_dx;
		ptr_velYdy = velY_dy;
		ptr_tmp = Tmp;
		ptr_lap_den = lap_den ;
*/
		ptr_StateVec = Umain;

	}if(grid == "MainGrid"){ // e vice-versa
/*		ptr_snd = vel_snd_arr_mid;
		ptr_den = den_mid;
		ptr_px = flxX_mid;
		ptr_py = flxY_mid;
		ptr_dendx = den_dx_mid;
		ptr_dendy = den_dy_mid;
		ptr_velXdx = velX_dx_mid;
		ptr_velXdy = velX_dy_mid;
		ptr_velYdx = velY_dx_mid;
		ptr_velYdy = velY_dy_mid;
		ptr_tmp = tmp_mid;
		ptr_lap_den = lap_den_mid ;
*/
		ptr_StateVec = Umid;
	}
}

float Fluid2D::SideAverage(const float * input_array, GridPoint2D p, char side){
	float avg;
	switch(side) {
		case 'N': avg=0.5f*(input_array[p.NE]+input_array[p.NW]);
			break;
		case 'S': avg=0.5f*(input_array[p.SE]+input_array[p.SW]);
			break;
		case 'E': avg=0.5f*(input_array[p.NE]+input_array[p.SE]);
			break;
		case 'W': avg=0.5f*(input_array[p.NW]+input_array[p.SW]);
			break;
		default: avg=0.0f;
	}
	return avg;
}



StateVec2D Fluid2D::SideAverage(const StateVec2D *input_array, GridPoint2D p, char side) {
	StateVec2D avg(input_array[p.C]);
	switch(side) {
		case 'N': avg=0.5f*(input_array[p.NE]+input_array[p.NW]);
			break;
		case 'S': avg=0.5f*(input_array[p.SE]+input_array[p.SW]);
			break;
		case 'E': avg=0.5f*(input_array[p.NE]+input_array[p.SE]);
			break;
		case 'W': avg=0.5f*(input_array[p.NW]+input_array[p.SW]);
			break;
		default: avg=input_array[p.C];
	}
	return avg;
}


float Fluid2D::DensitySource(StateVec2D U) {
	return 0;
}

float Fluid2D::XMomentumSource(StateVec2D U) {
	return 0;
}

float Fluid2D::YMomentumSource(StateVec2D U) {
	return 0;
}

float Fluid2D::TemperatureSource(StateVec2D U) {
	return 0;
}



float Fluid2D::TemperatureFluxX(StateVec2D U) {
	return 0;
}

float Fluid2D::TemperatureFluxY(StateVec2D U) {
	return 0;
}

void Fluid2D::CopyFields() {
	float mass;
	for (int i = 0; i < Nx*Ny; ++i) {
		Den[i]=Umain[i].n();
		mass= DensityToMass(Den[i]);
		VelX[i]=Umain[i].px()/mass;
		VelY[i]=Umain[i].py()/mass;
		Tmp[i] =Umain[i].tmp();
	}
}

