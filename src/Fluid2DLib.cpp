/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
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
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx * Ny]();
	velX_dx	= new float[Nx * Ny]();
	velX_dy	= new float[Nx * Ny]();
	velY_dx	= new float[Nx * Ny]();
	velY_dy	= new float[Nx * Ny]();

	den_dx = new float[Nx * Ny]();
	den_dy = new float[Nx * Ny]();

	lap_den = new float[Nx*Ny]();
	lap_den_mid =new float[(Nx-1)*(Ny-1)]();

	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve
    lap_tmp = new float[Nx*Ny]();

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	tmp_mid		= new float[(Nx-1)*(Ny-1)]();
	den_mid		= new float[(Nx-1)*(Ny-1)]();
	velX_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_mid	= new float[(Nx-1)*(Ny-1)]();
	velX_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velX_dy_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dx_mid	= new float[(Nx-1)*(Ny-1)]();
	velY_dy_mid	= new float[(Nx-1)*(Ny-1)]();

	den_dx_mid = new float[(Nx-1)*(Ny-1)]();
	den_dy_mid = new float[(Nx-1)*(Ny-1)]();

	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();
	vel_snd_arr_mid	= new float[(Nx-1)*(Ny-1)]();
}

Fluid2D::~Fluid2D() = default;

void Fluid2D::SetSound(){
	for(int kp=0; kp<=Nx*Ny-1; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		vel_snd_arr[kp]= Sound_Velocity_Anisotropy(i*dx, j*dy , vel_snd);
	}
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) { //correr todos os pontos da grelha secundaria
		div_t divresult;
		divresult = div(ks, Nx - 1);
		auto j = static_cast<float>(divresult.quot);
		auto i = static_cast<float>(divresult.rem);
		vel_snd_arr_mid[ks]= Sound_Velocity_Anisotropy((i+0.5f)*dx, (j+0.5f)*dy , vel_snd);
	}
}




void Fluid2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Den[c] = 1.0f + 0.005f * (noise - 0.5f);

		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        Tmp[c] =  .2f + 0.005f * (noise - 0.5f);
	}
}

void Fluid2D::InitialCondWave() {
	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
			//Den[i + j * Nx] = 1.0f+0.3f*sin(2.0f*MAT_PI*i*dx/lengX) ;
			Den[i + j * Nx] = 1.0f+0.05f*sin(5.0f*2.0f*MAT_PI*i*dx/lengX) ;
			//Den[i + j * Nx] = 1.0f+0.3f/cosh((i*dx-0.5f)*20.0f);
			VelX[i + j * Nx] = 0.0f;
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
			Den[i + j * Nx] = 1.0f + densi;
			VelX[i + j * Nx] = 0.1f;
		}
	}
}


void Fluid2D::MassFluxToVelocity(const string grid) {
	if(grid=="MainGrid"){
		for(int c=0; c <= Nx * Ny - 1; c++){
			VelX[c]= FlxX[c] / Den[c];
			VelY[c]= FlxY[c] / Den[c];
		}
	}
	if(grid=="MidGrid"){
		for(int c=0; c <= (Nx-1) * (Ny-1) - 1; c++){
			velX_mid[c]= flxX_mid[c] / den_mid[c];
			velY_mid[c]= flxY_mid[c] / den_mid[c];
		}
	}
}

void Fluid2D::VelocityToCurrent() {
	for(int c=0; c <= Nx * Ny - 1; c++){
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::Richtmyer(){
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			GridPoint midpoint(ks, Nx, Ny, true);

			float den_avg   = 0.25f * (Den[midpoint.SW] + Den[midpoint.SE] + Den[midpoint.NW] + Den[midpoint.NE]);
			float flx_x_avg = 0.25f * (FlxX[midpoint.SW] + FlxX[midpoint.SE] + FlxX[midpoint.NW] + FlxX[midpoint.NE]);
			float flx_y_avg = 0.25f * (FlxY[midpoint.SW] + FlxY[midpoint.SE] + FlxY[midpoint.NW] + FlxY[midpoint.NE]);
            float tmp_avg   = 0.25f * (Tmp[midpoint.SW] + Tmp[midpoint.SE] + Tmp[midpoint.NW] + Tmp[midpoint.NE]);
			float n_dx=0.0f,n_dy=0.0f;
			if(therm_diff) {
				n_dx =0.25f * (den_dx[midpoint.SW] + den_dx[midpoint.SE] + den_dx[midpoint.NW] + den_dx[midpoint.NE]);
				n_dy =0.25f * (den_dy[midpoint.SW] + den_dy[midpoint.SE] + den_dy[midpoint.NW] + den_dy[midpoint.NE]) ;
			}
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
			if(therm_diff){
				tmp_mid[ks] = tmp_avg
				              -0.5f * (dt / dx) * (TemperatureFluxX(midpoint, 'E') - TemperatureFluxX(midpoint, 'W'))
				              -0.5f * (dt / dy) * (TemperatureFluxY(midpoint, 'N') - TemperatureFluxY(midpoint, 'S'))
							  +0.5f * dt * TemperatureSource(den_avg, flx_x_avg, flx_y_avg, n_dx, n_dy, 0.0f, 0.0f);
			}
		}
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
	ChooseGridPointers("MainGrid");
#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			GridPoint mainpoint(kp, Nx, Ny, false);
			if( kp%Nx!=Nx-1 && kp%Nx!=0){
				float den_old = Den[kp];
				float flx_x_old = FlxX[kp];
				float flx_y_old = FlxY[kp];
                float tmp_old = Tmp[kp];

                Den[kp] = den_old - (dt/dx)*(DensityFluxX(mainpoint, 'E') - DensityFluxX(mainpoint, 'W'))
						          - (dt/dy)*(DensityFluxY(mainpoint, 'N') - DensityFluxY(mainpoint, 'S'))
						          + dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxX[kp] = flx_x_old - (dt/dx)*(XMomentumFluxX(mainpoint, 'E') - XMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(XMomentumFluxY(mainpoint, 'N') - XMomentumFluxY(mainpoint, 'S'))
						             + dt*XMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxY[kp] = flx_y_old - (dt/dx)*(YMomentumFluxX(mainpoint, 'E') - YMomentumFluxX(mainpoint, 'W'))
						             - (dt/dy)*(YMomentumFluxY(mainpoint, 'N') - YMomentumFluxY(mainpoint, 'S'))
				                     + dt*YMomentumSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				if(therm_diff) {
					Tmp[kp] = tmp_old - (dt / dx) * (TemperatureFluxX(mainpoint, 'E') - TemperatureFluxX(mainpoint, 'W'))
					          - (dt / dy) * (TemperatureFluxY(mainpoint, 'N') - TemperatureFluxY(mainpoint, 'S'))
							  + dt * TemperatureSource(den_old, flx_x_old, flx_y_old, den_dx[kp], den_dy[kp], 0.0f, 0.0f);
				}
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
		if(!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(FlxX[pos_end]) || !isfinite(FlxX[pos_ini])){
			cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
			CloseHdf5File();
			exit(EXIT_FAILURE);
		}
	data_preview << t << "\t"
	<< Den[pos_end]  << "\t"
	<< FlxX[pos_end] << "\t"
	<< Den[pos_ini]  << "\t"
	<< FlxX[pos_ini] << "\n";
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
		GridPoint point(kp,Nx,Ny,false);
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			lap_flxX[kp] = Laplacian19( point, VelX, kin_vis);
			lap_flxY[kp] = Laplacian19( point, VelY, kin_vis);
		}
	}
}

void Fluid2D::TemperatureLaplacianWeighted19() {
#pragma omp parallel for default(none) shared(lap_tmp,Tmp)
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
	    GridPoint point(kp,Nx,Ny,false);
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


void Fluid2D::SaveSound() {
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

	this->MassFluxToVelocity("MainGrid");
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
#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			flx_x_old = FlxX[kp];
			flx_y_old = FlxY[kp];
            tmp_old = Tmp[kp];

            FlxX[kp] = flx_x_old + lap_flxX[kp];
			FlxY[kp] = flx_y_old + lap_flxY[kp];
            Tmp[kp] = tmp_old + lap_tmp[kp];
        }
	}
}

/*
void Fluid2D::VelocityGradient() {
	int stride = Nx;
	float m_east,m_west,m_north,m_south;
#pragma omp parallel for default(none) private(m_east,m_west,m_north,m_south) shared(Nx,Ny,dx,dy,stride,FlxX,FlxY,Den,velX_dx,velX_dy,velY_dx,velY_dy)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){
		if( kp%stride!=stride-1 && kp%stride!=0){
			GridPoint point(kp,Nx,Ny,false);
			m_east =  DensityToMass(Den[point.E]);
			m_west =  DensityToMass(Den[point.W]);
			m_north = DensityToMass(Den[point.N]);
			m_south = DensityToMass(Den[point.S]);
			velX_dx[kp] = ( (FlxX[point.E]/m_east) -(FlxX[point.W]/m_west) )/(2.0f*dx);
			velX_dy[kp] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy);
			velY_dx[kp] = ( (FlxY[point.E]/m_east) -(FlxY[point.W]/m_west) )/(2.0f*dx);
			velY_dy[kp] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy);
		}
	}

	for(int i=1 ; i<=Nx-2; i++){
		float m_top,m_southsouth;
		int top= i + (Ny - 1) * stride;
		GridPoint point(top,Nx,Ny,false);
		int southsouth= i + (Ny - 3) * stride;
		m_top =   DensityToMass(Den[point.C]);
		m_east =  DensityToMass(Den[point.E]);
		m_west =  DensityToMass(Den[point.W]);
		m_south = DensityToMass(Den[point.S]);
		m_southsouth = DensityToMass(Den[southsouth]);
		velX_dx[top] = ( (FlxX[point.E]/m_east)-(FlxX[point.W]/m_west) )/(2.0f*dx);
		velX_dy[top] = ( 3.0f*(FlxX[point.C]/m_top) -4.0f*(FlxX[point.S]/m_south) +1.0f*(FlxX[southsouth]/m_southsouth) )/(2.0f*dy); //backward finite difference
		velY_dx[top] = ( (FlxY[point.E]/m_east)-(FlxY[point.W]/m_west) )/(2.0f*dx);
		velY_dy[top] = ( 3.0f*(FlxY[point.C]/m_top) -4.0f*(FlxY[point.S]/m_south) +1.0f*(FlxY[southsouth]/m_southsouth)  )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=Nx-2; i++){
		float m_bottom,m_northnorth;
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,false);
		int northnorth=i+2*stride;
		m_bottom = DensityToMass(Den[point.C]);
		m_east =   DensityToMass(Den[point.E]);
		m_west =   DensityToMass(Den[point.W]);
		m_north =  DensityToMass(Den[point.N]);
		m_northnorth = DensityToMass(Den[northnorth]);
		velX_dx[bottom] = ( (FlxX[point.E]/m_east)-(FlxX[point.W]/m_west) )/(2.0f*dx);
		velX_dy[bottom] = ( -3.0f*(FlxX[point.C]/m_bottom) +4.0f*(FlxX[point.N]/m_north)-1.0f*(FlxX[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
		velY_dx[bottom] = ( (FlxY[point.E]/m_east)-(FlxY[point.W]/m_west) )/(2.0f*dx);
		velY_dy[bottom] = ( -3.0f*(FlxY[point.C]/m_bottom) +4.0f*(FlxY[point.N]/m_north)-1.0f*(FlxY[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=Ny-2;j++){
		float m_left, m_easteast;
		int left = 0 + j*stride;
		GridPoint point(left,Nx,Ny,false);
		int easteast = left + 2;
		m_left =  DensityToMass(Den[point.C]);
		m_north = DensityToMass(Den[point.N]);
		m_south = DensityToMass(Den[point.S]);
		m_east =  DensityToMass(Den[point.E]);
		m_easteast = DensityToMass(Den[easteast]);
		velX_dx[left] = ( -3.0f*(FlxX[point.C]/m_left) +4.0f*(FlxX[point.E]/m_east)-1.0f*(FlxX[easteast]/m_easteast) )/(2.0f*dx); //forward difference
		velX_dy[left] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy);
		velY_dx[left] = ( -3.0f*(FlxY[point.C]/m_left) +4.0f*(FlxY[point.E]/m_east)-1.0f*(FlxY[easteast]/m_easteast) )/(2.0f*dx); //forward difference
		velY_dy[left] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy);
	}
	for(int j=1; j<=Ny-2;j++){
		float m_rigth, m_westwest;
		int right = (Nx-1) + j*stride;
		GridPoint point(right,Nx,Ny,false);
		int westwest = right-2;
		m_rigth = DensityToMass(Den[point.C]);
		m_north = DensityToMass(Den[point.N]);
		m_south = DensityToMass(Den[point.S]);
		m_west =  DensityToMass(Den[point.W]);
		m_westwest = DensityToMass(Den[westwest]);
		velX_dx[right] = ( 3.0f*(FlxX[point.C]/m_rigth) -4.0f*(FlxX[point.W]/m_west) +1.0f*(FlxX[westwest]/m_westwest) )/(2.0f*dx); //backwar difference
		velX_dy[right] = ( (FlxX[point.N]/m_north)-(FlxX[point.S]/m_south) )/(2.0f*dy);
		velY_dx[right] = ( 3.0f*(FlxY[point.C]/m_rigth) -4.0f*(FlxY[point.W]/m_west) +1.0f*(FlxY[westwest]/m_westwest) )/(2.0f*dx);
		velY_dy[right] = ( (FlxY[point.N]/m_north)-(FlxY[point.S]/m_south) )/(2.0f*dy);
	}
	int kp;
	// i=0 j=0 forward x forward y
	kp = 0 + 0*Nx;
	float m_0 =  DensityToMass(Den[kp]);
	float m_x1 = DensityToMass(Den[kp+1]);
	float m_x2 = DensityToMass(Den[kp+2]);
	float m_y1 = DensityToMass(Den[kp+Nx]);
	float m_y2 = DensityToMass(Den[kp+2*Nx]);
	velX_dx[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1]/m_x1 )-1.0f*(FlxX[kp+2]/m_x2 )  )/(2.0f*dx);
	velX_dy[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1*Nx]/m_y1 )-1.0f*(FlxX[kp+2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1]/m_x1 )-1.0f*(FlxY[kp+2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1*Nx]/m_y1 )-1.0f*(FlxY[kp+2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=(Nx-1) j=0 backward x forward y
	kp = (Nx-1) + 0*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp-1]);
	m_x2 = DensityToMass(Den[kp-2]);
	m_y1 = DensityToMass(Den[kp+Nx]);
	m_y2 = DensityToMass(Den[kp+2*Nx]);
	velX_dx[kp] = (  3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-1]/m_x1 )+1.0f*(FlxX[kp-2]/m_x2 )  )/(2.0f*dx);
	velX_dy[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+Nx]/m_y1 )-1.0f*(FlxX[kp+2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = (  3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-1]/m_x1 )+1.0f*(FlxY[kp-2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+Nx]/m_y1 )-1.0f*(FlxY[kp+2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=0 j=(Ny-1) forward x backward y
	kp = 0 + (Ny-1)*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp+1]);
	m_x2 = DensityToMass(Den[kp+2]);
	m_y1 = DensityToMass(Den[kp-Nx]);
	m_y2 = DensityToMass(Den[kp-2*Nx]);
	velX_dx[kp] = ( -3.0f*(FlxX[kp]/m_0 ) +4.0f*(FlxX[kp+1]/m_x1 )-1.0f*(FlxX[kp+2]/m_x1 )  )/(2.0f*dx);
	velX_dy[kp] = (  3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-Nx]/m_y1 )+1.0f*(FlxX[kp-2*Nx]/m_y2 )  )/(2.0f*dy);
	velY_dx[kp] = ( -3.0f*(FlxY[kp]/m_0 ) +4.0f*(FlxY[kp+1]/m_x1 )-1.0f*(FlxY[kp+2]/m_x2 )  )/(2.0f*dx);
	velY_dy[kp] = (  3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-Nx]/m_y1 )+1.0f*(FlxY[kp-2*Nx]/m_y2 )  )/(2.0f*dy);

	// i=(Nx-1) j=(Ny-1) backward x backward y
	kp = (Nx-1) + (Ny-1)*Nx;
	m_0 =  DensityToMass(Den[kp]);
	m_x1 = DensityToMass(Den[kp-1]);
	m_x2 = DensityToMass(Den[kp-2]);
	m_y1 = DensityToMass(Den[kp-Nx]);
	m_y2 = DensityToMass(Den[kp-2*Nx]);
	velX_dx[kp] = ( 3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-1]/m_x1 )+1.0f*(FlxX[kp-2]/m_x2 ) )/(2.0f*dx);
	velX_dy[kp] = ( 3.0f*(FlxX[kp]/m_0 ) -4.0f*(FlxX[kp-Nx]/m_y1 )+1.0f*(FlxX[kp-2*Nx]/m_y2 ) )/(2.0f*dy);
	velY_dx[kp] = ( 3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-1]/m_x1 )+1.0f*(FlxY[kp-2]/m_x2 ) )/(2.0f*dx);
	velY_dy[kp] = ( 3.0f*(FlxY[kp]/m_0 ) -4.0f*(FlxY[kp-Nx]/m_y1 )+1.0f*(FlxY[kp-2*Nx]/m_y2 ) )/(2.0f*dy);
}
*/

/*
void Fluid2D::VelocityGradientMid() {
	float m_east,m_west,m_north,m_south;
#pragma omp parallel for default(none) private(m_east,m_west,m_north,m_south) shared(Nx,Ny,dx,dy,flxX_mid,flxY_mid,den_mid,velX_dx_mid,velX_dy_mid,velY_dx_mid,velY_dy_mid)
	for(int ks=1+(Nx-1); ks<=(Nx-3)+(Ny-3)*(Nx-1); ks++){
		if( ks%(Nx-1)!=Nx-2 && ks%(Nx-1)!=0) {
			GridPoint point(ks,Nx,Ny,true);
			m_east = DensityToMass(den_mid[point.E]);
			m_west = DensityToMass(den_mid[point.W]);
			m_north = DensityToMass(den_mid[point.N]);
			m_south = DensityToMass(den_mid[point.S]);
			velX_dx_mid[ks] = ( (flxX_mid[point.E]/m_east) - (flxX_mid[point.W]/m_west) ) / (2.0f * dx);
			velX_dy_mid[ks] = ( (flxX_mid[point.N]/m_north) - (flxX_mid[point.S]/m_south) ) / (2.0f * dy);
			velY_dx_mid[ks] = ( (flxY_mid[point.E]/m_east) - (flxY_mid[point.W]/m_west) ) / (2.0f * dx);
			velY_dy_mid[ks] = ( (flxY_mid[point.N]/m_north) - (flxY_mid[point.S]/m_south) ) / (2.0f * dy);
		}
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){
		float m_top,m_southsouth;
		int top= i + ((Ny-1) - 1) * (Nx-1);
		GridPoint point(top,Nx,Ny,true);
		int southsouth= i + ((Ny-1) - 3) * (Nx-1);
		m_top = DensityToMass(den_mid[point.C]);
		m_east = DensityToMass(den_mid[point.E]);
		m_west = DensityToMass(den_mid[point.W]);
		m_south = DensityToMass(den_mid[point.S]);
		m_southsouth = DensityToMass(den_mid[southsouth]);
		velX_dx_mid[top] = ( (flxX_mid[point.E]/m_east)-(flxX_mid[point.W]/m_west) )/(2.0f*dx);
		velX_dy_mid[top] = ( 3.0f*(flxX_mid[point.C]/m_top) -4.0f*(flxX_mid[point.S]/m_south) +1.0f*(flxX_mid[southsouth]/m_southsouth) )/(2.0f*dy); //backward finite difference
		velY_dx_mid[top] = ( (flxY_mid[point.E]/m_east)-(flxY_mid[point.W]/m_west) )/(2.0f*dx);
		velY_dy_mid[top] = ( 3.0f*(flxY_mid[point.C]/m_top) -4.0f*(flxY_mid[point.S]/m_south) +1.0f*(flxY_mid[southsouth]/m_southsouth)  )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){ // fundo rede principal, ou seja j=0
		float m_bottom,m_northnorth;
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,true);
		int northnorth=i+2*(Nx-1);
		m_bottom = DensityToMass(den_mid[point.C]);
		m_east = DensityToMass(den_mid[point.E]);
		m_west = DensityToMass(den_mid[point.W]);
		m_north = DensityToMass(den_mid[point.N]);
		m_northnorth = DensityToMass(den_mid[northnorth]);
		velX_dx_mid[bottom] = ( (flxX_mid[point.E]/m_east)-(flxX_mid[point.W]/m_west) )/(2.0f*dx);
		velX_dy_mid[bottom] = ( -3.0f*(flxX_mid[point.C]/m_bottom) +4.0f*(flxX_mid[point.N]/m_north)-1.0f*(flxX_mid[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
		velY_dx_mid[bottom] = ( (flxY_mid[point.E]/m_east)-(flxY_mid[point.W]/m_west) )/(2.0f*dx);
		velY_dy_mid[bottom] = ( -3.0f*(flxY_mid[point.C]/m_bottom) +4.0f*(flxY_mid[point.N]/m_north)-1.0f*(flxY_mid[northnorth]/m_northnorth) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado esquerdo da rede principal ou seja i=0
		float m_left, m_easteast;
		int left = 0 + j*(Nx-1);
		GridPoint point(left,Nx,Ny,true);
		int easteast = left + 2;
		m_north = DensityToMass(den_mid[point.N]);
		m_south = DensityToMass(den_mid[point.S]);
		m_east = DensityToMass(den_mid[point.E]);
		m_easteast = DensityToMass(den_mid[easteast]);
		m_left= DensityToMass(den_mid[point.C]);
		velX_dx_mid[left] = ( -3.0f*(flxX_mid[point.C]/m_left) +4.0f*(flxX_mid[point.E]/m_east)-1.0f*(flxX_mid[easteast]/m_easteast)  )/(2.0f*dx); //forward difference
		velX_dy_mid[left] = ( (flxX_mid[point.N]/m_north)-(flxX_mid[point.S]/m_south) )/(2.0f*dy);
		velY_dx_mid[left] = ( -3.0f*(flxY_mid[point.C]/m_left) +4.0f*(flxY_mid[point.E]/m_east)-1.0f*(flxY_mid[easteast]/m_easteast) )/(2.0f*dx); //forward difference
		velY_dy_mid[left] = ( (flxY_mid[point.N]/m_north)-(flxY_mid[point.S]/m_south) )/(2.0f*dy);
	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado direito da rede principal ou seja i=((Nx-1)-1)
		float m_rigth, m_westwest;
		int right = ((Nx-1)-1) + j*(Nx-1);
		GridPoint point(right,Nx,Ny,true);
		int westwest = right-2;
		m_north = DensityToMass(den_mid[point.N]);
		m_south = DensityToMass(den_mid[point.S]);
		m_rigth= DensityToMass(den_mid[point.C]);
		m_west = DensityToMass(den_mid[point.W]);
		m_westwest = DensityToMass(den_mid[westwest]);
		velX_dx_mid[right] = ( 3.0f*(flxX_mid[point.C]/m_rigth) -4.0f*(flxX_mid[point.W]/m_west) +1.0f*(flxX_mid[westwest]/m_westwest) )/(2.0f*dx); //backwar difference
		velX_dy_mid[right] = ( (flxX_mid[point.N]/m_north)-(flxX_mid[point.S]/m_south) )/(2.0f*dy);
		velY_dx_mid[right] = ( 3.0f*(flxY_mid[point.C]/m_rigth) -4.0f*(flxY_mid[point.W]/m_west) +1.0f*(flxY_mid[westwest]/m_westwest) )/(2.0f*dx);
		velY_dy_mid[right] = ( (flxY_mid[point.N]/m_north)-(flxY_mid[point.S]/m_south) )/(2.0f*dy);
	}

	int ks;
// i=0 j=0 forward x forward y
	ks = 0 + 0*(Nx-1);
	float m_0 = DensityToMass(den_mid[ks]);
	float m_x1 = DensityToMass(den_mid[ks+1]);
	float m_x2 = DensityToMass(den_mid[ks+2]);
	float m_y1 = DensityToMass(den_mid[ks+(Nx-1)]);
	float m_y2 = DensityToMass(den_mid[ks+2*(Nx-1)]);
	velX_dx_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1]/m_x1 )-1.0f*(flxX_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velX_dy_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1*(Nx-1)]/m_y1 )-1.0f*(flxX_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1]/m_x1 )-1.0f*(flxY_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1*(Nx-1)]/m_y1 )-1.0f*(flxY_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=((Nx-1)-1) j=0 backward x forward y
	ks = ((Nx-1)-1) + 0*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks-1]);
	m_x2 = DensityToMass(den_mid[ks-2]);
	m_y1 = DensityToMass(den_mid[ks+(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks+2*(Nx-1)]);
	velX_dx_mid[ks] = (  3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-1]/m_x1 )+1.0f*(flxX_mid[ks-2]/m_x2 )  )/(2.0f*dx);
	velX_dy_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+(Nx-1)]/m_y1 )-1.0f*(flxX_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = (  3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-1]/m_x1 )+1.0f*(flxY_mid[ks-2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+(Nx-1)]/m_y1 )-1.0f*(flxY_mid[ks+2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=0 j=((Ny-1)-1) forward x backward y
	ks = 0 + ((Ny-1)-1)*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks+1]);
	m_x2 = DensityToMass(den_mid[ks+2]);
	m_y1 = DensityToMass(den_mid[ks-(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks-2*(Nx-1)]);
	velX_dx_mid[ks] = ( -3.0f*(flxX_mid[ks]/m_0 ) +4.0f*(flxX_mid[ks+1]/m_x1 )-1.0f*(flxX_mid[ks+2]/m_x1 )  )/(2.0f*dx);
	velX_dy_mid[ks] = (  3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxX_mid[ks-2*(Nx-1)]/m_y2 )  )/(2.0f*dy);
	velY_dx_mid[ks] = ( -3.0f*(flxY_mid[ks]/m_0 ) +4.0f*(flxY_mid[ks+1]/m_x1 )-1.0f*(flxY_mid[ks+2]/m_x2 )  )/(2.0f*dx);
	velY_dy_mid[ks] = (  3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxY_mid[ks-2*(Nx-1)]/m_y2 )  )/(2.0f*dy);

// i=((Nx-1)-1) j=((Ny-1)-1) backward x backward y
	ks = ((Nx-1)-1) + ((Ny-1)-1)*(Nx-1);
	m_0 = DensityToMass(den_mid[ks]);
	m_x1 = DensityToMass(den_mid[ks-1]);
	m_x2 = DensityToMass(den_mid[ks-2]);
	m_y1 = DensityToMass(den_mid[ks-(Nx-1)]);
	m_y2 = DensityToMass(den_mid[ks-2*(Nx-1)]);
	velX_dx_mid[ks] = ( 3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-1]/m_x1 )+1.0f*(flxX_mid[ks-2]/m_x2 ) )/(2.0f*dx);
	velX_dy_mid[ks] = ( 3.0f*(flxX_mid[ks]/m_0 ) -4.0f*(flxX_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxX_mid[ks-2*(Nx-1)]/m_y2 ) )/(2.0f*dy);
	velY_dx_mid[ks] = ( 3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-1]/m_x1 )+1.0f*(flxY_mid[ks-2]/m_x2 ) )/(2.0f*dx);
	velY_dy_mid[ks] = ( 3.0f*(flxY_mid[ks]/m_0 ) -4.0f*(flxY_mid[ks-(Nx-1)]/m_y1 )+1.0f*(flxY_mid[ks-2*(Nx-1)]/m_y2 ) )/(2.0f*dy);
}
*/

/*
void Fluid2D::DensityGradient() {
	int stride = Nx;
#pragma omp parallel for default(none) shared(Nx,Ny,dx,dy,stride,Den,den_dx,den_dy)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){
		if( kp%stride!=stride-1 && kp%stride!=0){
			GridPoint point(kp,Nx,Ny,false);
			den_dx[kp] = ( Den[point.E] - Den[point.W] )/(2.0f*dx);
			den_dy[kp] = ( Den[point.N] - Den[point.S] )/(2.0f*dy);
		}
	}

	for(int i=1 ; i<=Nx-2; i++){ // topo rede principal, ou seja j=(Ny - 1)
		int top= i + (Ny - 1) * stride;
		GridPoint point(top,Nx,Ny,false);
		int southsouth= i + (Ny - 3) * stride;
		den_dx[top] = ( (Den[point.E])-(Den[point.W]) )/(2.0f*dx); //OK
		den_dy[top] = ( 3.0f*(Den[point.C]) -4.0f*(Den[point.S]) +1.0f*(Den[southsouth]) )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=Nx-2; i++){ // fundo rede principal, ou seja j=0
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,false);
		int northnorth=i+2*stride;
		den_dx[bottom] = ( (Den[point.E])-(Den[point.W]) )/(2.0f*dx);
		den_dy[bottom] = ( -3.0f*(Den[point.C]) +4.0f*(Den[point.N])-1.0f*(Den[northnorth]) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=Ny-2;j++){ //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j*stride;
		GridPoint point(left,Nx,Ny,false);
		int easteast = left + 2;
		den_dx[left] = ( -3.0f*(Den[point.C]) +4.0f*(Den[point.E])-1.0f*(Den[easteast]) )/(2.0f*dx); //forward difference
		den_dy[left] = ( (Den[point.N])-(Den[point.S]) )/(2.0f*dy); //OK
	}
	for(int j=1; j<=Ny-2;j++){ //lado direito da rede principal ou seja i=(Nx-1)
		int right = (Nx-1) + j*stride;
		GridPoint point(right,Nx,Ny,false);
		int westwest = right-2;
		den_dx[right] = ( 3.0f*(Den[point.C]) -4.0f*(Den[point.W]) +1.0f*(Den[westwest]) )/(2.0f*dx); //backwar difference
		den_dy[right] = ( (Den[point.N])-(Den[point.S]) )/(2.0f*dy);
	}

	int kp;
	// i=0 j=0 forward x forward y
	kp = 0 + 0*Nx;
	den_dx[kp] = ( -3.0f*(Den[kp] ) +4.0f*(Den[kp+1] )-1.0f*(Den[kp+2] )  )/(2.0f*dx);
	den_dy[kp] = ( -3.0f*(Den[kp] ) +4.0f*(Den[kp+1*Nx] )-1.0f*(Den[kp+2*Nx] )  )/(2.0f*dy);
	// i=(Nx-1) j=0 backward x forward y
	kp = (Nx-1) + 0*Nx;
	den_dx[kp] = (  3.0f*(Den[kp] ) -4.0f*(Den[kp-1] )+1.0f*(Den[kp-2] )  )/(2.0f*dx);
	den_dy[kp] = ( -3.0f*(Den[kp] ) +4.0f*(Den[kp+Nx] )-1.0f*(Den[kp+2*Nx] )  )/(2.0f*dy);
	// i=0 j=(Ny-1) forward x backward y
	kp = 0 + (Ny-1)*Nx;
	den_dx[kp] = ( -3.0f*(Den[kp] ) +4.0f*(Den[kp+1] )-1.0f*(Den[kp+2] )  )/(2.0f*dx);
	den_dy[kp] = (  3.0f*(Den[kp] ) -4.0f*(Den[kp-Nx] )+1.0f*(Den[kp-2*Nx] )  )/(2.0f*dy);
	// i=(Nx-1) j=(Ny-1) backward x backward y
	kp = (Nx-1) + (Ny-1)*Nx;
	den_dx[kp] = ( 3.0f*(Den[kp] ) -4.0f*(Den[kp-1] )+1.0f*(Den[kp-2] ) )/(2.0f*dx);
	den_dy[kp] = ( 3.0f*(Den[kp] ) -4.0f*(Den[kp-Nx] )+1.0f*(Den[kp-2*Nx] ) )/(2.0f*dy);
}
*/

/*
void Fluid2D::DensityGradientMid() {

#pragma omp parallel for default(none)  shared(Nx,Ny,dx,dy,den_mid,den_dx_mid,den_dy_mid)
	for(int ks=1+(Nx-1); ks<=(Nx-3)+(Ny-3)*(Nx-1); ks++){ //correr todos os pontos da grelha secundaria de _mid EVITANDO FRONTEIRAS
		if( ks%(Nx-1)!=Nx-2 && ks%(Nx-1)!=0) {
			GridPoint point(ks,Nx,Ny,true);
			den_dx_mid[ks] = ( (den_mid[point.E]) - (den_mid[point.W]) ) / (2.0f * dx);
			den_dy_mid[ks] = ( (den_mid[point.N]) - (den_mid[point.S]) ) / (2.0f * dy);
		}
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){ // topo rede principal, ou seja j=((Ny-1) - 1)
		int top= i + ((Ny-1) - 1) * (Nx-1);
		GridPoint point(top,Nx,Ny,true);
		int southsouth= i + ((Ny-1) - 3) * (Nx-1);
		den_dx_mid[top] = ( (den_mid[point.E])-(den_mid[point.W]) )/(2.0f*dx);
		den_dy_mid[top] = ( 3.0f*(den_mid[point.C]) -4.0f*(den_mid[point.S]) +1.0f*(den_mid[southsouth]) )/(2.0f*dy); //backward finite difference
	}
	for(int i=1 ; i<=(Nx-1)-2; i++){ // fundo rede principal, ou seja j=0
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,true);
		int northnorth=i+2*(Nx-1);
		den_dx_mid[bottom] = ( (den_mid[point.E])-(den_mid[point.W]) )/(2.0f*dx);
		den_dy_mid[bottom] = ( -3.0f*(den_mid[point.C]) +4.0f*(den_mid[point.N])-1.0f*(den_mid[northnorth]) )/(2.0f*dy); //forward finite difference
	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j*(Nx-1);
		GridPoint point(left,Nx,Ny,true);
		int easteast = left + 2;
		den_dx_mid[left] = ( -3.0f*(den_mid[point.C]) +4.0f*(den_mid[point.E])-1.0f*(den_mid[easteast])  )/(2.0f*dx); //forward difference
		den_dy_mid[left] = ( (den_mid[point.N])-(den_mid[point.S]) )/(2.0f*dy);
	}
	for(int j=1; j<=(Ny-1)-2;j++){ //lado direito da rede principal ou seja i=((Nx-1)-1)
		int right = ((Nx-1)-1) + j*(Nx-1);
		GridPoint point(right,Nx,Ny,true);
		int westwest = right-2;
		den_dx_mid[right] = ( 3.0f*(den_mid[point.C]) -4.0f*(den_mid[point.W]) +1.0f*(den_mid[westwest]) )/(2.0f*dx); //backwar difference
		den_dy_mid[right] = ( (den_mid[point.N])-(den_mid[point.S]) )/(2.0f*dy);
	}
//os 4 cantos em que ambas a derivadas nao podem ser centradas
	int ks;
// i=0 j=0 forward x forward y
	ks = 0 + 0*(Nx-1);
	den_dx_mid[ks] = ( -3.0f*(den_mid[ks] ) +4.0f*(den_mid[ks+1] )-1.0f*(den_mid[ks+2] )  )/(2.0f*dx);
	den_dy_mid[ks] = ( -3.0f*(den_mid[ks] ) +4.0f*(den_mid[ks+1*(Nx-1)] )-1.0f*(den_mid[ks+2*(Nx-1)] )  )/(2.0f*dy);
// i=((Nx-1)-1) j=0 backward x forward y
	ks = ((Nx-1)-1) + 0*(Nx-1);
	den_dx_mid[ks] = (  3.0f*(den_mid[ks] ) -4.0f*(den_mid[ks-1] )+1.0f*(den_mid[ks-2] )  )/(2.0f*dx);
	den_dy_mid[ks] = ( -3.0f*(den_mid[ks] ) +4.0f*(den_mid[ks+(Nx-1)] )-1.0f*(den_mid[ks+2*(Nx-1)] )  )/(2.0f*dy);
// i=0 j=((Ny-1)-1) forward x backward y
	ks = 0 + ((Ny-1)-1)*(Nx-1);
	den_dx_mid[ks] = ( -3.0f*(den_mid[ks] ) +4.0f*(den_mid[ks+1] )-1.0f*(den_mid[ks+2] )  )/(2.0f*dx);
	den_dy_mid[ks] = (  3.0f*(den_mid[ks] ) -4.0f*(den_mid[ks-(Nx-1)] )+1.0f*(den_mid[ks-2*(Nx-1)] )  )/(2.0f*dy);
// i=((Nx-1)-1) j=((Ny-1)-1) backward x backward y
	ks = ((Nx-1)-1) + ((Ny-1)-1)*(Nx-1);
	den_dx_mid[ks] = ( 3.0f*(den_mid[ks] ) -4.0f*(den_mid[ks-1] )+1.0f*(den_mid[ks-2] ) )/(2.0f*dx);
	den_dy_mid[ks] = ( 3.0f*(den_mid[ks] ) -4.0f*(den_mid[ks-(Nx-1)] )+1.0f*(den_mid[ks-2*(Nx-1)] ) )/(2.0f*dy);
}
*/

/*
void Fluid2D::DensityLaplacian(){
	int stride = Nx;
#pragma omp parallel for default(none) shared(Nx,Ny,dx,dy,stride,Den,lap_den)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){
		if( kp%stride!=stride-1 && kp%stride!=0){
			GridPoint point(kp,Nx,Ny,false);
			lap_den[kp] = (-4.0f*Den[point.C] +  Den[point.N] + Den[point.S] + Den[point.E] + Den[point.W] )/(dx*dx);
		}
	}

	// -2	5	-4	1
	for(int i=1 ; i<=Nx-2; i++){ // topo rede principal, ou seja j=(Ny - 1)
		int top= i + (Ny - 1) * stride;
		GridPoint point(top,Nx,Ny,false);
		int southsouth= i + (Ny - 3) * stride;
		int southsouthsouth= i + (Ny - 4) * stride;
		float aux1 = -2.0f*Den[point.C] + Den[point.E] + Den[point.W]; //OK
		float aux2 = -2.0f*Den[point.C] +5.0f*Den[point.S] -4.0f*Den[southsouth] +1.0f*Den[southsouthsouth];
		lap_den[top] = (aux1+aux2)/(dx*dx);
	}

	// 2	−5	4	−1
	for(int i=1 ; i<=Nx-2; i++){ // fundo rede principal, ou seja j=0
		int bottom=i; //i+0*nx
		GridPoint point(bottom,Nx,Ny,false);
		int northnorth=i+2*stride;
		int northnorthnorth=i+3*stride;
		float aux1 = -2.0f*Den[point.C] + Den[point.E] + Den[point.W]; //OK
		float aux2 =  2.0f*Den[point.C] -5.0f*Den[point.S] +4.0f*Den[northnorth] -1.0f*Den[northnorthnorth];
		lap_den[bottom] = (aux1+aux2)/(dx*dx);
	}

	// 2	−5	4	−1
	for(int j=1; j<=Ny-2;j++){ //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j*stride;
		GridPoint point(left,Nx,Ny,false);
		int easteast = left + 2;
		int easteasteast = left + 3;
		float aux1 =  -2.0f*Den[point.C] + Den[point.N] + Den[point.S];
		float aux2 =  2.0f*Den[point.C] -5.0f*Den[point.E] +4.0f*Den[easteast] -1.0f*Den[easteasteast];
		lap_den[left] = (aux1+aux2)/(dx*dx);
	}

	// -2	5	-4	1
	for(int j=1; j<=Ny-2;j++){ //lado direito da rede principal ou seja i=(Nx-1)
		int right = (Nx-1) + j*stride;
		GridPoint point(right,Nx,Ny,false);
		int westwest = right-2;
		int westwestwest = right-3;
		float aux1 =  -2.0f*Den[point.C] + Den[point.N] + Den[point.S];
		float aux2 =  -2.0f*Den[point.C] +5.0f*Den[point.W] -4.0f*Den[westwest] +1.0f*Den[westwestwest];
		lap_den[right] = (aux1+aux2)/(dx*dx);
	}


	//2	−5	4	−1
	int kp;
	// i=0 j=0 forward x forward y
	kp = 0 + 0*Nx;
	float aux1 = 2.0f*Den[kp] -5.0f*Den[kp+1] +4.0f*Den[kp+2] -1.0f*Den[kp+3];
	float aux2 = 2.0f*Den[kp] -5.0f*Den[kp+Nx] +4.0f*Den[kp+2*Nx] -1.0f*Den[kp+3*Nx];
	lap_den[kp] = (aux1+aux2)/(dx*dx);
	// i=(Nx-1) j=0 backward x forward y
	kp = (Nx-1) + 0*Nx;
	aux1 = -2.0f*Den[kp] +5.0f*Den[kp-1] -4.0f*Den[kp-2] +1.0f*Den[kp-3];
	aux2 = 2.0f*Den[kp] -5.0f*Den[kp+Nx] +4.0f*Den[kp+2*Nx] -1.0f*Den[kp+3*Nx];
	lap_den[kp] =(aux1+aux2)/(dx*dx);
	// i=0 j=(Ny-1) forward x backward y
	kp = 0 + (Ny-1)*Nx;
	aux1 = 2.0f*Den[kp] -5.0f*Den[kp+1] +4.0f*Den[kp+2] -1.0f*Den[kp+3];
	aux2 = -2.0f*Den[kp] +5.0f*Den[kp-Nx] -4.0f*Den[kp-2*Nx] +1.0f*Den[kp-3*Nx];
	lap_den[kp] = (aux1+aux2)/(dx*dx);
	// i=(Nx-1) j=(Ny-1) backward x backward y
	kp = (Nx-1) + (Ny-1)*Nx;
	aux1 = -2.0f*Den[kp] +5.0f*Den[kp-1] -4.0f*Den[kp-2] +1.0f*Den[kp-3];
	aux2 = -2.0f*Den[kp] +5.0f*Den[kp-Nx] -4.0f*Den[kp-2*Nx] +1.0f*Den[kp-3*Nx];
	lap_den[kp] = (aux1+aux2)/(dx*dx);

}
*/



float Fluid2D::DensityToMass(float density) {
	return density;
}

float Fluid2D::DensityFluxX(GridPoint p, char side) {
//	float * px_ptr;
	float px;
/*	if(p.IsMidGrid){
		px_ptr = FlxX;
	}else{
		px_ptr = flxX_mid;
	}
	if (side == 'E'){
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
	}
	if (side == 'W'){
		px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
	}*/
	px = SideAverage(ptr_px,p,side);
	return px;
}

float Fluid2D::DensityFluxY(GridPoint p, char side) {
	//float * py_ptr;
	float py;
	/*
	if(p.IsMidGrid){
		py_ptr = FlxY;
	}else{
		py_ptr = flxY_mid;
	}
	if (side == 'N'){
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
	}
	if (side == 'S'){
		py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
	}*/
	py = SideAverage(ptr_py,p,side);
	return py;
}

float Fluid2D::XMomentumFluxX(GridPoint p, char side) {
//	float * den_ptr;
//	float * px_ptr;
	float den;
	float px;
/*
	if(p.IsMidGrid){
		den_ptr = Den;
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
 */
	den = SideAverage(ptr_den,p,side);
	px = SideAverage(ptr_px,p,side);
	return px * px / den + den;
}

float Fluid2D::XMomentumFluxY(GridPoint p, char side) {
	//float * den_ptr;
	//float * px_ptr;
	//float * py_ptr;

	float den;
	float px;
	float py;
/*
	if(p.IsMidGrid){
		den_ptr = Den;
		px_ptr = FlxX;
		py_ptr = FlxY;
	}else{
		den_ptr = den_mid;
		px_ptr = flxX_mid;
		py_ptr = flxY_mid;
	}
	if (side == 'N'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.NW]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.NW]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
	}
	if (side == 'S'){
		den = 0.5f*(den_ptr[p.SE] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.SE] + px_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
	}*/
	den = SideAverage(ptr_den,p,side);
	px = SideAverage(ptr_px,p,side);
	py = SideAverage(ptr_py,p,side);
	return px * py / den;
}


float Fluid2D::YMomentumFluxY(GridPoint p, char side) {
	//float * den_ptr;
	//float * py_ptr;
	float den;
	float py;
	den = SideAverage(ptr_den,p,side);

	py = SideAverage(ptr_py,p,side);
	/*
	if(p.IsMidGrid){
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
	}*/
	return py * py / den + den;
}

float Fluid2D::YMomentumFluxX(GridPoint p, char side) {
	//float * den_ptr;
	//float * px_ptr;
	//float * py_ptr;
	float den;
	float px;
	float py;
	den = SideAverage(ptr_den,p,side);
	px = SideAverage(ptr_px,p,side);
	py = SideAverage(ptr_py,p,side);
	/*
	if(p.IsMidGrid){
		den_ptr = Den;
		px_ptr = FlxX;
		py_ptr = FlxY;
	}else{
		den_ptr = den_mid;
		px_ptr = flxX_mid;
		py_ptr = flxY_mid;
	}
	if (side == 'E'){
		den = 0.5f*(den_ptr[p.NE] + den_ptr[p.SE]);
		px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
		py = 0.5f*(py_ptr[p.NE] + py_ptr[p.SE]);
	}
	if (side == 'W'){
		den = 0.5f*(den_ptr[p.NW] + den_ptr[p.SW]);
		px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
		py = 0.5f*(py_ptr[p.NW] + py_ptr[p.SW]);
	}*/
	return px * py / den;
}


float Fluid2D::TemperatureFluxX(GridPoint p, char side) {
    //float * px_ptr;
	//float * den_ptr;
	//float * tmp_ptr;
	float tmp;
    float px;
	float den;

	den = SideAverage(ptr_den,p,side);
	px = SideAverage(ptr_px,p,side);
	tmp = SideAverage(ptr_tmp,p,side);
	/*
    if(p.IsMidGrid){
	    den_ptr = Den;
        px_ptr = FlxX;
	    tmp_ptr = Tmp;
    }else{
	    den_ptr = den_mid;
        px_ptr = flxX_mid;
	    tmp_ptr = tmp_mid;
    }
    if (side == 'E'){
	    den = 0.5f*(den_ptr[p.NE] + den_ptr[p.SE]);
        px = 0.5f*(px_ptr[p.NE] + px_ptr[p.SE]);
	    tmp = 0.5f*(tmp_ptr[p.NE] + tmp_ptr[p.SE]);
    }
    if (side == 'W'){
	    den = 0.5f*(den_ptr[p.NW] + den_ptr[p.SW]);
        px = 0.5f*(px_ptr[p.NW] + px_ptr[p.SW]);
	    tmp = 0.5f*(tmp_ptr[p.NW] + tmp_ptr[p.SW]);
    }*/


    //return px * vel_fer * vel_fer + px / DensityToMass(den);
	return PHYS_FERMI_CNVC*px  + tmp*px / DensityToMass(den);
}


float Fluid2D::TemperatureFluxY(GridPoint p, char side) {
    //float * py_ptr;
    float py;
	//float * tmp_ptr;
	float tmp;
	//float * den_ptr;
	float den;

	den = SideAverage(ptr_den,p,side);
	py = SideAverage(ptr_py,p,side);
	tmp = SideAverage(ptr_tmp,p,side);


	/*
    if(p.IsMidGrid){
	    den_ptr = Den;
        py_ptr = FlxY;
        tmp_ptr=Tmp;
    }else{
	    den_ptr = den_mid;
        py_ptr = flxY_mid;
        tmp_ptr=tmp_mid;
    }
    if (side == 'N'){
	    den = 0.5f*(den_ptr[p.NE] + den_ptr[p.NW]);
        py = 0.5f*(py_ptr[p.NE] + py_ptr[p.NW]);
	    tmp = 0.5f*(tmp_ptr[p.NE] + tmp_ptr[p.NW]);
    }
    if (side == 'S'){
	    den = 0.5f*(den_ptr[p.SE] + den_ptr[p.SW]);
        py = 0.5f*(py_ptr[p.SE] + py_ptr[p.SW]);
	    tmp = 0.5f*(tmp_ptr[p.SE] + tmp_ptr[p.SW]);
    }
    */

    //return 0.0f*py * vel_fer * vel_fer  + py / DensityToMass(den) ;
	return PHYS_FERMI_CNVC*py  + tmp*py / DensityToMass(den) ;
}



float Fluid2D::Laplacian19(GridPoint p, float *input_ptr, float constant) {
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

void Fluid2D::ChooseGridPointers(string grid) {
	if(grid == "MidGrid"){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		ptr_snd = vel_snd_arr;
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
	}if(grid == "MainGrid"){ // e vice-versa
		ptr_snd = vel_snd_arr_mid;
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
	}
}

float Fluid2D::SideAverage(const float * input_array, GridPoint p,char side){
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
