// 2D version



#include "includes/Fluid2DLib.h"
#include "SetUpParametersLib.h"
#include "GrapheneFluid2DLib.h"


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

	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2fwc=%.2f", vel_snd, vel_fer, kin_vis, col_freq,cyc_freq);
	file_infix = buffer;
	// main grid variables Nx*Ny
	Den 		= new float[Nx * Ny]();
	VelX 		= new float[Nx * Ny]();
	VelY 		= new float[Nx * Ny]();
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx * Ny]();


	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();  
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
//#pragma omp parallel for default(none) shared(Nx,Ny,maxrand,rd)
	for (int c = 0; c < Nx*Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		Den[c] = 1.0f + 0.005f * (noise - 0.5f);
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


void Fluid2D::MassFluxToVelocity(){
//#pragma omp parallel for default(none) shared(VelX,VelY,FlxX,FlxY,Den)
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] / Den[c];
		VelY[c]= FlxY[c] / Den[c];
		//CurX[c] = VelX[c] * Den[c];
		//CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::VelocityToCurrent() {
//#pragma omp parallel for default(none) shared(Nx,Ny,VelX,VelY,CurX,CurY,Den)
	for(int c=0; c <= Nx * Ny - 1; c++){
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::Richtmyer(){


#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr,dt,dx)
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
			int northeast,northwest,southeast,southwest;
			float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
			float  sound_north, sound_south ,sound_east ,sound_west;

			div_t divresult;
			divresult = div (ks,Nx-1);
			int j=divresult.quot;
			int i=divresult.rem;

			northeast=i+1+(j+1)*Nx;
			northwest=i+(j+1)*Nx;
			southeast=i+1+j*Nx;
			southwest=i+j*Nx;

			sound_north = 0.5f * (vel_snd_arr[northeast] + vel_snd_arr[northwest]);
			sound_south = 0.5f*(vel_snd_arr[southeast] + vel_snd_arr[southwest]);
			sound_east = 0.5f*(vel_snd_arr[northeast] + vel_snd_arr[southeast]);
			sound_west = 0.5f*(vel_snd_arr[northwest] + vel_snd_arr[southwest]);

			den_north = 0.5f*(Den[northeast] + Den[northwest]);
			den_south = 0.5f*(Den[southeast] + Den[southwest]);
			den_east = 0.5f*(Den[northeast] + Den[southeast]);
			den_west = 0.5f*(Den[northwest] + Den[southwest]);

			px_north = 0.5f*(FlxX[northeast] + FlxX[northwest]);
			px_south = 0.5f*(FlxX[southeast] + FlxX[southwest]);
			px_east = 0.5f*(FlxX[northeast] + FlxX[southeast]);
			px_west = 0.5f*(FlxX[northwest] + FlxX[southwest]);

			py_north = 0.5f*(FlxY[northeast] + FlxY[northwest]);
			py_south = 0.5f*(FlxY[southeast] + FlxY[southwest]);
			py_east = 0.5f*(FlxY[northeast] + FlxY[southeast]);
			py_west = 0.5f*(FlxY[northwest] + FlxY[southwest]);

			m_east=sqrt(den_east*den_east*den_east);
			m_west=sqrt(den_west*den_west*den_west);
			m_north=sqrt(den_north*den_north*den_north);
			m_south=sqrt(den_south*den_south*den_south);

			float den_avg = 0.25f * (Den[southwest] + Den[southeast] + Den[northwest] + Den[northeast]);
			float flx_x_avg = 0.25f * (FlxX[southwest] + FlxX[southeast] + FlxX[northwest] + FlxX[northeast]);
			float flx_y_avg = 0.25f * (FlxY[southwest] + FlxY[southeast] + FlxY[northwest] + FlxY[northeast]);

			den_mid[ks] = den_avg
					-0.5f*(dt/dx)*(
						DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
						DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
						DensityFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*DensitySource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);


			flxX_mid[ks] = flx_x_avg
					-0.5f*(dt/dx)*(
						MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
						MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
						MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*MassFluxXSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);

			flxY_mid[ks] = flx_y_avg
					-0.5f*(dt/dx)*(
						MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
						MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
					-0.5f*(dt/dy)*(
						MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
						MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))//;
					+0.5f*dt*MassFluxYSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
		}

#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr_mid,dt,dx)
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
			int northeast,northwest,southeast,southwest;
			float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
			float  sound_north, sound_south ,sound_east ,sound_west;

			div_t divresult;
			divresult = div (kp,Nx);
			int j=divresult.quot;
			int i=divresult.rem;


			if( kp%Nx!=Nx-1 && kp%Nx!=0){
				northeast=i+j*(Nx-1);
				northwest=i-1+j*(Nx-1);
				southeast=i+(j-1)*(Nx-1);
				southwest=i-1+(j-1)*(Nx-1);

				sound_north = 0.5f * (vel_snd_arr_mid[northeast] + vel_snd_arr_mid[northwest]);
				sound_south = 0.5f*(vel_snd_arr_mid[southeast] + vel_snd_arr_mid[southwest]);
				sound_east = 0.5f*(vel_snd_arr_mid[northeast] + vel_snd_arr_mid[southeast]);
				sound_west = 0.5f*(vel_snd_arr_mid[northwest] + vel_snd_arr_mid[southwest]);

				den_north = 0.5f*(den_mid[northeast]+den_mid[northwest]);
				den_south = 0.5f*(den_mid[southeast]+den_mid[southwest]);
				den_east = 0.5f*(den_mid[northeast]+den_mid[southeast]);
				den_west = 0.5f*(den_mid[northwest]+den_mid[southwest]);

				px_north = 0.5f*(flxX_mid[northeast]+flxX_mid[northwest]);
				px_south = 0.5f*(flxX_mid[southeast]+flxX_mid[southwest]);
				px_east = 0.5f*(flxX_mid[northeast]+flxX_mid[southeast]);
				px_west = 0.5f*(flxX_mid[northwest]+flxX_mid[southwest]);
				
				py_north = 0.5f*(flxY_mid[northeast]+flxY_mid[northwest]);
				py_south = 0.5f*(flxY_mid[southeast]+flxY_mid[southwest]);
				py_east = 0.5f*(flxY_mid[northeast]+flxY_mid[southeast]);
				py_west = 0.5f*(flxY_mid[northwest]+flxY_mid[southwest]);

				m_east=sqrt(den_east*den_east*den_east);
				m_west=sqrt(den_west*den_west*den_west);
				m_north=sqrt(den_north*den_north*den_north);
				m_south=sqrt(den_south*den_south*den_south);

				float den_old = Den[kp];
				float flx_x_old = FlxX[kp];
				float flx_y_old = FlxY[kp];


				Den[kp] = den_old - (dt / dx) * (
							DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
							DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
							DensityFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxX[kp] = flx_x_old - (dt / dx) * (
							MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
							MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
							MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*MassFluxXSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
				FlxY[kp] = flx_y_old - (dt / dx) * (
							MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
							MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
						-(dt/dy)*(
							MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
							MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))
						+dt*MassFluxYSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);

			}
		}
}

void Fluid2D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0f;
}



float  Fluid2D::DensityFluxX(__attribute__((unused)) float n, float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s){
	float f_1;
	f_1 = flx_x;
	return f_1;
}
float  Fluid2D::DensityFluxY(__attribute__((unused)) float n, __attribute__((unused)) float flx_x, float flx_y, __attribute__((unused)) float mass, __attribute__((unused)) float s){
	float f_1;
	f_1 = flx_y;
	return f_1;
}

float  Fluid2D::MassFluxXFluxX(float n,float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_2;
	f_2 = flx_x * flx_x / n + n;
	return f_2;
}
float  Fluid2D::MassFluxXFluxY(float n,float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_2;
	f_2 = flx_x * flx_y / n;
	return f_2;
}
float  Fluid2D::MassFluxYFluxX(float n,float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_3;
	f_3 = flx_x * flx_y / n;
	return f_3;
}
float  Fluid2D::MassFluxYFluxY(float n, __attribute__((unused)) float flx_x, float flx_y,__attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_3;
	f_3 = flx_y * flx_y / n + n;
	return f_3;
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
	this->MassFluxToVelocity();

//calculate laplacians
#pragma omp parallel for  default(none) shared(lap_flxX,lap_flxY,VelX,VelY,Nx,Ny)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
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
	this->MassFluxToVelocity();
	float sx=kin_vis*dt/(dx*dx);
	float sy=kin_vis*dt/(dy*dy);
#pragma omp parallel for default(none) shared(lap_flxX,lap_flxY,VelX,VelY,sx,sy)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		int north, south, east, west, northeast, northwest,southeast, southwest ;
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
			northeast = i + 1 + (j + 1) * Nx;
			northwest = i - 1 + (j + 1) * Nx;
			southeast = i + 1 + (j - 1) * Nx;
			southwest = i - 1 + (j - 1) * Nx;
			lap_flxX[kp] = (4.0f*sx*sy-2.0f*sx-2.0f*sy)*VelX[kp] +
			               sx*sy*( VelX[northeast] + VelX[southeast] + VelX[northwest] + VelX[southwest])
			               + sy*(1.0f-2.0f*sx)*(VelX[north] + VelX[south])
			               + sx*(1.0f-2.0f*sy)*(VelX[west] + VelX[east]);
			lap_flxY[kp] = (4.0f*sx*sy-2.0f*sx-2.0f*sy)*VelY[kp] +
			               sx*sy*( VelY[northeast] + VelY[southeast] + VelY[northwest] + VelY[southwest])
			               + sy*(1.0f-2.0f*sx)*(VelY[north] + VelY[south])
			               + sx*(1.0f-2.0f*sy)*(VelY[west] + VelY[east]);
		}
	}
}


void Fluid2D:: ParabolicOperatorFtcs() {
	VelocityLaplacianFtcs();
	ForwardTimeOperator();
}

void Fluid2D::ParabolicOperatorWeightedExplicit19() {
	VelocityLaplacianWeighted19();
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

	this->MassFluxToVelocity();
	this->VelocityToCurrent();
	string str_time = to_string(TimeStepCounter / snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

	DataSet dataset_den = GrpDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", HDF5FLOAT, atr_dataspace);
	float currenttime=static_cast<float>(TimeStepCounter) * dt;
	atr_step_den.write(HDF5INT, &TimeStepCounter);
	atr_time_den.write(HDF5FLOAT , &currenttime);
	atr_step_den.close();
	atr_time_den.close();

	dataset_den.write(Den, HDF5FLOAT);
	dataset_den.close();

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

float Fluid2D::MassFluxXSource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}
float Fluid2D::MassFluxYSource(__attribute__((unused)) float n,__attribute__((unused)) float flx_x,__attribute__((unused)) float flx_y,__attribute__((unused)) float mass,__attribute__((unused)) float s) {
	return 0.0f;
}

void Fluid2D::ForwardTimeOperator() {
#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,lap_flxX,lap_flxY,dt,cyc_freq)
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		float flx_x_old,flx_y_old;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			flx_x_old=FlxX[kp];
			flx_y_old=FlxY[kp];
			FlxX[kp] = flx_x_old + lap_flxX[kp];
			FlxY[kp] = flx_y_old + lap_flxY[kp];
		}
	}
}
/*


void Fluid2D::RichtmyerFirstStep(){
#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr,dt,dx)
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
		int northeast,northwest,southeast,southwest;
		float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
		float  sound_north, sound_south ,sound_east ,sound_west;

		div_t divresult;
		divresult = div (ks,Nx-1);
		int j=divresult.quot;
		int i=divresult.rem;

		northeast=i+1+(j+1)*Nx;
		northwest=i+(j+1)*Nx;
		southeast=i+1+j*Nx;
		southwest=i+j*Nx;

		sound_north = 0.5f * (vel_snd_arr[northeast] + vel_snd_arr[northwest]);
		sound_south = 0.5f*(vel_snd_arr[southeast] + vel_snd_arr[southwest]);
		sound_east = 0.5f*(vel_snd_arr[northeast] + vel_snd_arr[southeast]);
		sound_west = 0.5f*(vel_snd_arr[northwest] + vel_snd_arr[southwest]);

		den_north = 0.5f*(Den[northeast] + Den[northwest]);
		den_south = 0.5f*(Den[southeast] + Den[southwest]);
		den_east = 0.5f*(Den[northeast] + Den[southeast]);
		den_west = 0.5f*(Den[northwest] + Den[southwest]);

		px_north = 0.5f*(FlxX[northeast] + FlxX[northwest]);
		px_south = 0.5f*(FlxX[southeast] + FlxX[southwest]);
		px_east = 0.5f*(FlxX[northeast] + FlxX[southeast]);
		px_west = 0.5f*(FlxX[northwest] + FlxX[southwest]);

		py_north = 0.5f*(FlxY[northeast] + FlxY[northwest]);
		py_south = 0.5f*(FlxY[southeast] + FlxY[southwest]);
		py_east = 0.5f*(FlxY[northeast] + FlxY[southeast]);
		py_west = 0.5f*(FlxY[northwest] + FlxY[southwest]);

		m_east=sqrt(den_east*den_east*den_east);
		m_west=sqrt(den_west*den_west*den_west);
		m_north=sqrt(den_north*den_north*den_north);
		m_south=sqrt(den_south*den_south*den_south);

		float den_avg = 0.25f * (Den[southwest] + Den[southeast] + Den[northwest] + Den[northeast]);
		float flx_x_avg = 0.25f * (FlxX[southwest] + FlxX[southeast] + FlxX[northwest] + FlxX[northeast]);
		float flx_y_avg = 0.25f * (FlxY[southwest] + FlxY[southeast] + FlxY[northwest] + FlxY[northeast]);

		den_mid[ks] = den_avg
		              -0.5f*(dt/dx)*(
				DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
				DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
		              -0.5f*(dt/dy)*(
				DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
				DensityFluxY(den_south, px_south, py_south,m_south,sound_south))//;
		              +0.5f*dt*DensitySource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);


		flxX_mid[ks] = flx_x_avg
		               -0.5f*(dt/dx)*(
				MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
				MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
		               -0.5f*(dt/dy)*(
				MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
				MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))//;
		               +0.5f*dt*MassFluxXSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);

		flxY_mid[ks] = flx_y_avg
		               -0.5f*(dt/dx)*(
				MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
				MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
		               -0.5f*(dt/dy)*(
				MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
				MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))//;
		               +0.5f*dt*MassFluxYSource(den_avg, flx_x_avg, flx_y_avg, 0.0f, 0.0f);
	}


}

void Fluid2D::RichtmyerSecondStep(){

#pragma omp parallel for default(none) shared(Nx,Ny,FlxX,FlxY,Den,flxX_mid,flxY_mid,den_mid,vel_snd_arr_mid,dt,dx)
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		int northeast,northwest,southeast,southwest;
		float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
		float  sound_north, sound_south ,sound_east ,sound_west;

		div_t divresult;
		divresult = div (kp,Nx);
		int j=divresult.quot;
		int i=divresult.rem;


		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			northeast=i+j*(Nx-1);
			northwest=i-1+j*(Nx-1);
			southeast=i+(j-1)*(Nx-1);
			southwest=i-1+(j-1)*(Nx-1);

			sound_north = 0.5f * (vel_snd_arr_mid[northeast] + vel_snd_arr_mid[northwest]);
			sound_south = 0.5f*(vel_snd_arr_mid[southeast] + vel_snd_arr_mid[southwest]);
			sound_east = 0.5f*(vel_snd_arr_mid[northeast] + vel_snd_arr_mid[southeast]);
			sound_west = 0.5f*(vel_snd_arr_mid[northwest] + vel_snd_arr_mid[southwest]);

			den_north = 0.5f*(den_mid[northeast]+den_mid[northwest]);
			den_south = 0.5f*(den_mid[southeast]+den_mid[southwest]);
			den_east = 0.5f*(den_mid[northeast]+den_mid[southeast]);
			den_west = 0.5f*(den_mid[northwest]+den_mid[southwest]);

			px_north = 0.5f*(flxX_mid[northeast]+flxX_mid[northwest]);
			px_south = 0.5f*(flxX_mid[southeast]+flxX_mid[southwest]);
			px_east = 0.5f*(flxX_mid[northeast]+flxX_mid[southeast]);
			px_west = 0.5f*(flxX_mid[northwest]+flxX_mid[southwest]);

			py_north = 0.5f*(flxY_mid[northeast]+flxY_mid[northwest]);
			py_south = 0.5f*(flxY_mid[southeast]+flxY_mid[southwest]);
			py_east = 0.5f*(flxY_mid[northeast]+flxY_mid[southeast]);
			py_west = 0.5f*(flxY_mid[northwest]+flxY_mid[southwest]);

			m_east=sqrt(den_east*den_east*den_east);
			m_west=sqrt(den_west*den_west*den_west);
			m_north=sqrt(den_north*den_north*den_north);
			m_south=sqrt(den_south*den_south*den_south);

			float den_old = Den[kp];
			float flx_x_old = FlxX[kp];
			float flx_y_old = FlxY[kp];
			Den[kp] = den_old - (dt / dx) * (
					DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
					DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
			          -(dt/dy)*(
					DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
					DensityFluxY(den_south, px_south, py_south,m_south,sound_south))
			          +dt*DensitySource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
			FlxX[kp] = flx_x_old - (dt / dx) * (
					MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
					MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
			           -(dt/dy)*(
					MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
					MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))
			           +dt*MassFluxXSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
			FlxY[kp] = flx_y_old - (dt / dx) * (
					MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
					MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
			           -(dt/dy)*(
					MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
					MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))
			           +dt*MassFluxYSource(den_old, flx_x_old, flx_y_old, 0.0f, 0.0f);
		}
	}
}

*/