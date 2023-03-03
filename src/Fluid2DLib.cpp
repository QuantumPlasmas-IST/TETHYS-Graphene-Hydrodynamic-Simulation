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
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();


	vel_snd_arr	= new float[Nx * Ny]();

	Umain = new StateVec2D[Nx * Ny]();
	Umid = new StateVec2D[(Nx - 1) * (Ny - 1)]();

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
		}
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




void Fluid2D::Richtmyer(){
	if(odd_vis){
		VelocityGradient(Umain,Nx,Ny);
	}
	RichtmyerStep1();
	if(odd_vis){
		VelocityGradient(Umid,Nx-1,Ny-1);
	}
	RichtmyerStep2();
}



void Fluid2D::RichtmyerStep1(){

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ChooseGridPointers("MidGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
		GridPoint2D midpoint(ks, Nx, Ny, true);

		StateVec2D Uavg(Umain[ks]);
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
		                -0.5f*(dt/dy)*(DensityFluxY(UNorth) - DensityFluxY(USouth))
						+0.5f*dt* DensitySource(Uavg);

		Umid[ks].px() = Uavg.px()
		                -0.5f*(dt/dx)*(XMomentumFluxX(UEast) - XMomentumFluxX(UWest))
		                -0.5f*(dt/dy)*(XMomentumFluxY(UNorth) - XMomentumFluxY(USouth))
		                +0.5f*dt*XMomentumSource(Uavg);

		Umid[ks].py() = Uavg.py()
		                -0.5f*(dt/dx)*(YMomentumFluxX(UEast) - YMomentumFluxX(UWest))
		                -0.5f*(dt/dy)*(YMomentumFluxY(UNorth) - YMomentumFluxY(USouth))
		                +0.5f*dt*YMomentumSource(Uavg);
	}



}


void Fluid2D::RichtmyerStep2(){

	ChooseGridPointers("MainGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		GridPoint2D mainpoint(kp, Nx, Ny, false);
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
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
			                //+ dt*EleDensitySource(Uold);
			Umain[kp].px() = Uold.px()
			                 - (dt/dx)*(XMomentumFluxX(UEast) - XMomentumFluxX(UWest))
			                 - (dt/dy)*(XMomentumFluxY(UNorth) - XMomentumFluxY(USouth));
			                 //+ dt*EleXMomentumSource(Uold);

			Umain[kp].py() = Uold.py()
			                 - (dt/dx)*(YMomentumFluxX(UEast) - YMomentumFluxX(UWest))
			                 - (dt/dy)*(YMomentumFluxY(UNorth) - YMomentumFluxY(USouth));
			                 //+ dt*EleYMomentumSource(Uold);
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
		if(!isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].px()) || !isfinite(Umain[pos_end].px())){
			cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
			CloseHdf5File();
			exit(EXIT_FAILURE);
		}
	data_preview << t <<"\t"<< Umain[pos_ini] <<"\t"<<Umain[pos_end]<< "\n";
}

void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}

void Fluid2D::VelocityLaplacianFtcs() {
	this->MassFluxToVelocity();
//#pragma omp parallel for  default(none) shared(lap_flxX,lap_flxY,VelX,VelY,Nx,Ny)
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
			Umain[kp].d2vx() =
					kin_vis*dt*(-4.0f * VelX[kp]  + VelX[north]  + VelX[south]  + VelX[east]  +
							VelX[west] ) / (dx * dx);
			Umain[kp].d2vy() =
					kin_vis*dt*(-4.0f * VelY[kp]  + VelY[north]  + VelY[south]  + VelY[east]  +
							VelY[west] ) / (dx * dx);
		}
	}
}

void Fluid2D::VelocityLaplacianWeighted19() {


	for (int kp = 0; kp < Nx * Ny ; kp++) {
		float den =Umain[kp].n();
		float mass = DensityToMass(den);
		VelX[kp] = Umain[kp].px()/mass;
		VelY[kp] = Umain[kp].py()/mass;
	}

#pragma omp parallel for default(none) shared(Umain,VelX,VelY)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		GridPoint2D point(kp,Nx,Ny,false);
		if (kp % Nx != Nx - 1 && kp % Nx != 0){
			//lap_flxX[kp] = Laplacian19( point, VelX, kin_vis);
			//lap_flxY[kp] = Laplacian19( point, VelY, kin_vis);
			Umain[kp].d2vx() = Laplacian19( point, VelX, kin_vis);
			Umain[kp].d2vy() = Laplacian19( point, VelY, kin_vis);
		}
	}
}

void Fluid2D::TemperatureLaplacianWeighted19() {
	for (int kp = 0; kp < Nx * Ny ; kp++) {
		Tmp[kp] = Umain[kp].tmp();
	}
#pragma omp parallel for default(none) shared(Umain,Tmp)
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
	    GridPoint2D point(kp,Nx,Ny,false);
        if (kp % Nx != Nx - 1 && kp % Nx != 0){
	        Umain[kp].d2tmp() = Laplacian19( point, Tmp, therm_diff );
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
//		vel_snd_arr[i]=vel_snd;
	}

	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocity", HDF5FLOAT, *DataspaceVelSnd);
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

	if(therm_diff!=0){
		DataSet dataset_tmp = GrpTmp->createDataSet(name_dataset, HDF5FLOAT, *DataspaceTmp);
		Attribute atr_step_tmp = dataset_tmp.createAttribute("time step", HDF5INT, atr_dataspace);
		Attribute atr_time_tmp = dataset_tmp.createAttribute("time", HDF5FLOAT, atr_dataspace);
		dataset_tmp.write(Tmp, HDF5FLOAT);
		dataset_tmp.close();
		atr_step_tmp.write(HDF5INT, &TimeStepCounter);
		atr_time_tmp.write(HDF5FLOAT, &currenttime);
		atr_step_tmp.close();
		atr_time_tmp.close();
	}
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



void Fluid2D::ForwardTimeOperator() {
//#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			flx_x_old = Umain[kp].px();
			flx_y_old = Umain[kp].py();
            tmp_old = Umain[kp].tmp();

			Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
			Umain[kp].py() = flx_y_old + Umain[kp].d2vy();
			Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();


		}
	}
}
void Fluid2D::ForwardTimeOperator(char field) {
//#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
		float flx_x_old, flx_y_old, tmp_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			switch(field) {
				case 'T': 	tmp_old = Umain[kp].tmp();
							Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();
							break;
				case 'V': 	flx_x_old = Umain[kp].px();
							flx_y_old = Umain[kp].py();
							Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
							Umain[kp].py() = flx_y_old + Umain[kp].d2vy();
							break;
				default: ;
			}
		}
	}
}

float Fluid2D::DensityToMass(float density) {
	return density;
}

float Fluid2D::DensityFluxX(StateVec2D U) {
	return U.px();
}



float Fluid2D::DensityFluxY(StateVec2D U) {
	return U.py();
}
float Fluid2D::XMomentumFluxX(StateVec2D U) {
	return U.px()*U.px()/U.n() + U.n();
}

float Fluid2D::XMomentumFluxY(StateVec2D U) {
	return U.px()*U.py()/U.n();
}


float Fluid2D::YMomentumFluxY(StateVec2D U) {
	return U.py()*U.py()/U.n() + U.n();
}


float Fluid2D::YMomentumFluxX(StateVec2D U) {
	return U.py()*U.px()/U.n();
}



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



void Fluid2D::ChooseGridPointers(const string &grid) {
	if(grid == "MidGrid"){
		ptr_StateVec = Umain;
	}if(grid == "MainGrid"){
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

void Fluid2D::VelocityGradient(StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
//#pragma omp parallel for default(none) shared(size_x, size_y, dx, dy, stride, array_in, array_out_x, array_out_y)
	for (int kp = 1 + size_x; kp <= size_x * size_y - size_x - 2; kp++) {
		if (kp % stride != stride - 1 && kp % stride != 0) {
			float mE = DensityToMass(Uarray[kp+1].n());
			float mW = DensityToMass(Uarray[kp-1].n());
			float mN = DensityToMass(Uarray[kp+stride].n());
			float mS = DensityToMass(Uarray[kp-stride].n());
			float vxE = Uarray[kp+1].px()/mE;
			float vxW = Uarray[kp-1].px()/mW;
			float vxN = Uarray[kp+stride].px()/mN;
			float vxS = Uarray[kp-stride].px()/mS;
			float vyE = Uarray[kp+1].py()/mE;
			float vyW = Uarray[kp-1].py()/mW;
			float vyN = Uarray[kp+stride].py()/mN;
			float vyS = Uarray[kp-stride].py()/mS;
			Uarray[kp].dxvx() = (vxE - vxW) / (2.0f * dx);
			Uarray[kp].dyvx() = (vxN - vxS) / (2.0f * dy);
			Uarray[kp].dxvy() = (vyE - vyW) / (2.0f * dx);
			Uarray[kp].dyvy() = (vyN - vyS) / (2.0f * dy);
		}
	}

	for (int i = 1; i <= size_x - 2; i++) { // topo rede principal, ou seja j=(size_y - 1)
		int top = i + (size_y - 1) * stride;
		int southsouth = i + (size_y - 3) * stride;
		float mC = DensityToMass(Uarray[top].n());
		float mE = DensityToMass(Uarray[top+1].n());
		float mW = DensityToMass(Uarray[top-1].n());
		float mS = DensityToMass(Uarray[top-stride].n());
		float mSS = DensityToMass(Uarray[southsouth].n());
		float vxC = Uarray[top].px()/mC;
		float vxE = Uarray[top+1].px()/mE;
		float vxW = Uarray[top-1].px()/mW;
		float vxSS = Uarray[southsouth].px()/mSS;
		float vxS = Uarray[top-stride].px()/mS;
		float vyC = Uarray[top].py()/mC;
		float vyE = Uarray[top+1].py()/mE;
		float vyW = Uarray[top-1].py()/mW;
		float vySS = Uarray[southsouth].py()/mSS;
		float vyS = Uarray[top-stride].py()/mS;
		Uarray[top].dxvx() = (vxE - vxW) / (2.0f * dx);
		Uarray[top].dxvy() = (vyE - vyW) / (2.0f * dx);
		Uarray[top].dyvx() = (3.0f * vxC - 4.0f * vxS + vxSS) /(2.0f * dy); //backward finite difference
		Uarray[top].dyvx() = (3.0f * vyC - 4.0f * vyS + vySS) /(2.0f * dy); //backward finite difference
	}

	for (int i = 1; i <= size_x - 2; i++) { // fundo rede principal, ou seja j=0
		int bottom = i; //i+0*nx
		int northnorth = i + 2 * stride;
		float mC = DensityToMass(Uarray[bottom].n());
		float mE = DensityToMass(Uarray[bottom+1].n());
		float mW = DensityToMass(Uarray[bottom-1].n());
		float mN = DensityToMass(Uarray[bottom+stride].n());
		float mNN = DensityToMass(Uarray[northnorth].n());
		float vxC = Uarray[bottom].px()/mC;
		float vxE = Uarray[bottom+1].px()/mE;
		float vxW = Uarray[bottom-1].px()/mW;
		float vxNN = Uarray[northnorth].px()/mNN;
		float vxN = Uarray[bottom+stride].px()/mN;
		float vyC = Uarray[bottom].py()/mC;
		float vyE = Uarray[bottom+1].py()/mE;
		float vyW = Uarray[bottom-1].py()/mW;
		float vyNN = Uarray[northnorth].py()/mNN;
		float vyN = Uarray[bottom+stride].py()/mN;

		Uarray[bottom].dxvx() = (vxE - vxW) / (2.0f * dx);
		Uarray[bottom].dxvy() = (vyE - vyW) / (2.0f * dx);
		Uarray[bottom].dyvx() = (-3.0f * vxC + 4.0f * vxN - vxNN) /(2.0f * dy); //backward finite difference
		Uarray[bottom].dyvx() = (-3.0f * vyC + 4.0f * vyN - vyNN) /(2.0f * dy); //backward finite difference
	}
	
	for (int j = 1; j <= size_y - 2; j++) { //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j * stride;
		int easteast = left + 2;


		float mC = DensityToMass(Uarray[left].n());
		float mE = DensityToMass(Uarray[left+1].n());
		float mEE = DensityToMass(Uarray[easteast].n());
		float mN = DensityToMass(Uarray[left+stride].n());
		float mS = DensityToMass(Uarray[left-stride].n());
		float vxC = Uarray[left].px()/mC;
		float vxE = Uarray[left+1].px()/mE;
		float vxEE = Uarray[easteast].px()/mEE;
		float vxN = Uarray[left+stride].px()/mN;
		float vxS = Uarray[left-stride].px()/mS;
		float vyC = Uarray[left].px()/mC;
		float vyE = Uarray[left+1].py()/mE;
		float vyEE = Uarray[easteast].py()/mEE;
		float vyN = Uarray[left+stride].py()/mN;
		float vyS = Uarray[left-stride].py()/mS;

		Uarray[left].dxvx() = (-3.0f * vxC+ 4.0f * vxE - vxEE) /(2.0f * dx); //forward difference
		Uarray[left].dxvy() = (-3.0f * vyC+ 4.0f * vyE - vyEE) /(2.0f * dx); //forward difference
		Uarray[left].dyvx() = (vxN - vxS) / (2.0f * dy); //OK
		Uarray[left].dyvy() = (vyN - vyS) / (2.0f * dy); //OK
	}

	for (int j = 1; j <= size_y - 2; j++) { //lado direito da rede principal ou seja i=(size_x-1)
		int right = (size_x - 1) + j * stride;
		int westwest = right - 2;


		float mC = DensityToMass(Uarray[right].n());
		float mWW = DensityToMass(Uarray[westwest].n());
		float mW = DensityToMass(Uarray[right-1].n());
		float mN = DensityToMass(Uarray[right+stride].n());
		float mS = DensityToMass(Uarray[right-stride].n());
		float vxC = Uarray[right].px()/mC;
		float vxWW = Uarray[westwest].px()/mWW;
		float vxW = Uarray[right-1].px()/mW;
		float vxN = Uarray[right+stride].px()/mN;
		float vxS = Uarray[right-stride].px()/mS;
		float vyC = Uarray[right].px()/mC;
		float vyWW = Uarray[westwest].py()/mWW;
		float vyW = Uarray[right-1].py()/mW;
		float vyN = Uarray[right+stride].py()/mN;
		float vyS = Uarray[right-stride].py()/mS;

		Uarray[right].dyvx() = (vxN - vxS) / (2.0f * dy); //OK
		Uarray[right].dyvy() = (vyN - vyS) / (2.0f * dy); //OK
		Uarray[right].dxvx() = (3.0f * vxC- 4.0f * vxW + vxWW) /(2.0f * dx); //backwar difference
		Uarray[right].dxvy() = (3.0f * vyC- 4.0f * vyW + vyWW) /(2.0f * dx);//backwar difference
	}

	int kp;
// i=0 j=0 forward x forward y
	kp = 0 + 0 * size_x;
	float mC = DensityToMass(Uarray[kp].n());
	float mE = DensityToMass(Uarray[kp+1].n());
	float mEE = DensityToMass(Uarray[kp+2].n());
	float mS = DensityToMass(Uarray[kp+stride].n());
	float mSS = DensityToMass(Uarray[kp+2*stride].n());
	float vxC = Uarray[kp].px()/mC;
	float vxE = Uarray[kp+1].px()/mE;
	float vxEE = Uarray[kp+2].px()/mEE;
	float vxS = Uarray[kp+stride].px()/mS;
	float vxSS = Uarray[kp+2*stride].px()/mSS;
	float vyC = Uarray[kp].px()/mC;
	float vyE = Uarray[kp+1].px()/mE;
	float vyEE = Uarray[kp+2].px()/mEE;
	float vyS = Uarray[kp+stride].px()/mS;
	float vySS = Uarray[kp+2*stride].px()/mSS;
	Uarray[kp].dxvx() = (-3.0f * vxC + 4.0f * vxE - vxEE ) / (2.0f * dx);
	Uarray[kp].dyvx() = (-3.0f * vxC + 4.0f * vxS - vxSS ) / (2.0f * dy);
	Uarray[kp].dxvy() = (-3.0f * vyC + 4.0f * vyE - vyEE ) / (2.0f * dx);
	Uarray[kp].dyvy() = (-3.0f * vyC + 4.0f * vyS - vySS ) / (2.0f * dy);
//-----------------------------------------------------------------------
// i=(size_x-1) j=0 backward x forward y
	kp = (size_x - 1) + 0 * size_x;
	mC = DensityToMass(Uarray[kp].n());
	float mW = DensityToMass(Uarray[kp-1].n());
	float mWW = DensityToMass(Uarray[kp-2].n());
	mS = DensityToMass(Uarray[kp+stride].n());
	mSS = DensityToMass(Uarray[kp+2*stride].n());
	vxC = Uarray[kp].px()/mC;
	float vxW = Uarray[kp-1].px()/mW;
	float vxWW = Uarray[kp-2].px()/mWW;
	vxS = Uarray[kp+stride].px()/mS;
	vxSS = Uarray[kp+2*stride].px()/mSS;
	vyC = Uarray[kp].px()/mC;
	float vyW = Uarray[kp-1].px()/mW;
	float vyWW = Uarray[kp-2].px()/mWW;
	vyS = Uarray[kp+stride].px()/mS;
	vySS = Uarray[kp+2*stride].px()/mSS;
	Uarray[kp].dxvx() = (3.0f * vxC - 4.0f * vxW + vxWW ) / (2.0f * dx);
	Uarray[kp].dyvx() = (-3.0f * vxC + 4.0f * vxS - vxSS ) / (2.0f * dy);
	Uarray[kp].dxvy() = (3.0f * vyC - 4.0f * vyW + vyWW ) / (2.0f * dx);
	Uarray[kp].dyvy() = (-3.0f * vyC + 4.0f * vyS - vySS ) / (2.0f * dy);
//-----------------------------------------------------------------------

// i=0 j=(size_y-1) forward x backward y
	kp = 0 + (size_y - 1) * size_x;
	mC = DensityToMass(Uarray[kp].n());
	mE = DensityToMass(Uarray[kp+1].n());
	mEE = DensityToMass(Uarray[kp+2].n());
	float mN = DensityToMass(Uarray[kp-stride].n());
	float mNN = DensityToMass(Uarray[kp-2*stride].n());
	vxC = Uarray[kp].px()/mC;
	vxE = Uarray[kp+1].px()/mE;
	vxEE = Uarray[kp+2].px()/mEE;
	float vxN = Uarray[kp-stride].px()/mN;
	float vxNN = Uarray[kp-2*stride].px()/mNN;
	vyC = Uarray[kp].px()/mC;
	vyE = Uarray[kp+1].px()/mE;
	vyEE = Uarray[kp+2].px()/mEE;
	float vyN = Uarray[kp-stride].px()/mN;
	float vyNN = Uarray[kp-2*stride].px()/mNN;
	Uarray[kp].dxvx() = (-3.0f * vxC + 4.0f * vxE - vxEE ) / (2.0f * dx);
	Uarray[kp].dyvx() = (3.0f * vxC - 4.0f * vxN + vxNN ) / (2.0f * dy);
	Uarray[kp].dxvy() = (-3.0f * vyC + 4.0f * vyE - vyEE ) / (2.0f * dx);
	Uarray[kp].dyvy() = (3.0f * vyC - 4.0f * vyN + vyNN ) / (2.0f * dy);

//-----------------------------------------------------------------------

// i=(size_x-1) j=(size_y-1) backward x backward y
	kp = (size_x - 1) + (size_y - 1) * size_x;
	mC = DensityToMass(Uarray[kp].n());
	mW = DensityToMass(Uarray[kp+1].n());
	mWW = DensityToMass(Uarray[kp+2].n());
	mN = DensityToMass(Uarray[kp-stride].n());
	mNN = DensityToMass(Uarray[kp-2*stride].n());
	vxC = Uarray[kp].px()/mC;
	vxW = Uarray[kp+1].px()/mW;
	vxWW = Uarray[kp+2].px()/mWW;
	vxN = Uarray[kp-stride].px()/mN;
	vxNN = Uarray[kp-2*stride].px()/mNN;
	vyC = Uarray[kp].px()/mC;
	vyW = Uarray[kp+1].px()/mW;
	vyWW = Uarray[kp+2].px()/mWW;
	vyN = Uarray[kp-stride].px()/mN;
	vyNN = Uarray[kp-2*stride].px()/mNN;
	Uarray[kp].dyvx() = (3.0f * vxC - 4.0f * vxN + vxNN ) / (2.0f * dy);
	Uarray[kp].dyvy() = (3.0f * vyC - 4.0f * vyN + vyNN ) / (2.0f * dy);
	Uarray[kp].dxvx() = (3.0f * vxC - 4.0f * vxW + vxWW ) / (2.0f * dx);
	Uarray[kp].dxvy() = (3.0f * vyC - 4.0f * vyW + vyWW ) / (2.0f * dx);

}

void Fluid2D::MassFluxToVelocity() {
	for (int i = 0; i < Nx*Ny; ++i) {
		float mass= DensityToMass(Umain[i].n());
		VelX[i]=Umain[i].px()/mass;
		VelY[i]=Umain[i].py()/mass;
	}
}

void Fluid2D::VelocityToCurrent() {
	for (int i = 0; i < Nx*Ny; ++i) {
		float mass= DensityToMass(Umain[i].n());
		CurX[i]=Den[i]*VelX[i];
		CurY[i]=Den[i]*VelY[i];
	}
}


