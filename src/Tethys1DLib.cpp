#include "Tethys1DLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905
#endif

using namespace H5;
using namespace std;


GrapheneFluid1D::GrapheneFluid1D(int size_n,  SetUpInput &input_parameters): Fluid1D(size_n, input_parameters){
	vel_fer = input_parameters.FermiVelocity;
	col_freq = input_parameters.CollisionFrequency;
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}

/*....................................................................*/	
/*.......... 1 Dimensional Fluid Class ...............................*/	
/*....................................................................*/	
Fluid1D::Fluid1D(int size_nx, const SetUpInput &input_parameters): TethysBase{size_nx, 0, 1}{
	Nx = size_nx;
	vel_snd = input_parameters.SoundVelocity;
	kin_vis = input_parameters.ShearViscosity;
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	Den = new float[size_nx]();
	Vel = new float[size_nx]();
	GradVel= new float[size_nx]();
	Cur = new float[size_nx]();
	DenCor = new float[size_nx]();
	VelCor = new float[size_nx]();
	CurCor = new float[size_nx]();
	den_mid = new float[size_nx - 1]();
	vel_mid = new float[size_nx - 1]();
	grad_vel_mid = new float[size_nx - 1]();
	vel_snd_arr = new float[size_nx - 1]();
}	
	
Fluid1D::~Fluid1D(){
	delete [] Den;
	delete [] Vel ;
	delete [] Cur ;
	delete [] den_mid ;
	delete [] vel_mid ;
	delete [] DenCor ;
	delete [] VelCor ;
	delete [] CurCor ;
	delete [] vel_snd_arr ;
	delete [] GradVel ;
	delete [] grad_vel_mid ;
}

float  Fluid1D::DensityFlux(float n,float v, __attribute__((unused)) float s){
	float f_1;
	f_1 = n * v;
	return f_1;
}
float  Fluid1D::VelocityFlux(float n,float v,float dv, __attribute__((unused)) float s){
	float f_2;
	f_2 = 0.5f * v * v + n - kin_vis * dv;
	return f_2;
}
float  Fluid1D::DensitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0;
}
float  Fluid1D::VelocitySource(__attribute__((unused)) float n,__attribute__((unused)) float v,__attribute__((unused)) float s){
	return 0;
}

void Fluid1D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dt = dx/10.0f;
}

void Fluid1D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx-1  ;i++){
		vel_snd_arr[i]= Sound_Velocity_Anisotropy(i*dx, vel_snd);
	}
}
		
void Fluid1D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) rd.max();

	for (int i = 0; i < Nx; i++ ){
		float noise = (float) rd()/ maxrand ;
		Den[i] = 1.0f + 0.005f * (noise - 0.5f);
	}
}

void Fluid1D::InitialCondTest(){
	//float mean=dx*Nx/2.0f;
	//float sigma=dx*Nx/10.0f;
	for (int i = 0; i < Nx; i++ ){
		//Vel[i] = 0.4f*exp(-0.5f*((i*dx-mean)*(i*dx-mean)/(sigma*sigma)))/sigma;
		Vel[i] = 1.0f+tanh(10.0f*(dx*i-0.5f));
	}
}

void Fluid1D::Smooth(int width){
	Average_Filter(Den, DenCor, Nx, width);
	Average_Filter(Vel, VelCor, Nx, width);
	Average_Filter(Cur, CurCor, Nx, width);
}



void Fluid1D::CreateFluidFile(){
	std::string previewfile = "preview_1D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid1D::WriteFluidFile(float t){
	try {
		int pos_end = Nx - 1 ;
		int pos_ini = 0;
		if (!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(Vel[pos_end]) ||
		    !isfinite(Vel[pos_ini])) {
			throw "ERROR: numerical method failed to converge";
		}
//data_preview << t << "\t" << DenCor[Nx - 1] << "\t" << VelCor[Nx - 1] << "\t" << DenCor[0] << "\t" << VelCor[0] << "\n";
		data_preview << t << "\t" << Den[pos_end] << "\t" << Vel[pos_end] << "\t" << Den[pos_ini] << "\t" << Vel[pos_ini] << "\n";
	}catch (const char* msg) {
		cerr << msg <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
}



void Fluid1D::Richtmyer(){
		//
		//Calculating the velocity gradient at k time
		//
		for ( int i = 1; i <= Nx-2 ; i++ )
		{
			GradVel[i] = (-0.5f * Vel[i - 1] + 0.5f * Vel[i + 1]) / dx;
		}
		GradVel[0] = (-1.5f * Vel[0] + 2.0f * Vel[1] - 0.5f * Vel[2]) / dx;
		GradVel[Nx - 1] = (0.5f * Vel[Nx - 1 - 2] - 2.0f * Vel[Nx - 1 - 1] + 1.5f * Vel[Nx - 1]) / dx;
		//
		//Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i <= Nx - 2; i++ )
		{
			den_mid[i] = 0.5f*(Den[i] + Den[i + 1] )
				- ( 0.5f*dt/dx ) * (DensityFlux(Den[i + 1], Vel[i + 1], vel_snd_arr[i]) - DensityFlux(Den[i], Vel[i], vel_snd_arr[i]) )
					+ ( 0.5f*dt    ) * DensitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
			vel_mid[i] = 0.5f*(Vel[i] + Vel[i + 1] )
				- ( 0.5f*dt/dx ) * (VelocityFlux(Den[i + 1], Vel[i + 1], GradVel[i + 1], vel_snd_arr[i]) - VelocityFlux(Den[i], Vel[i], GradVel[i], vel_snd_arr[i]) )
					+ ( 0.5f*dt    ) * VelocitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
		}
		//
		//  Calculating the velocity gradient at k+1/2 time
		//
		for ( int i = 1; i <= Nx-3 ; i++ )
		{
			grad_vel_mid[i] =(-0.5f*vel_mid[i-1]+0.5f*vel_mid[i+1])/dx;
		}
		grad_vel_mid[0] = (-1.5f*vel_mid[0]+2.0f*vel_mid[1]-0.5f*vel_mid[2])/dx;
		grad_vel_mid[(Nx-1)-1] = ( 0.5f*vel_mid[(Nx-1)-3]-2.0f*vel_mid[(Nx-1)-2]+1.5f*vel_mid[(Nx-1)-1])/dx;
		//
		// Remaining step 
		//
		for ( int i = 1; i <= Nx - 2; i++ )
		{
			float den_old = Den[i];
			float vel_old = Vel[i];
			Den[i] = Den[i] - (dt / dx) * (DensityFlux(den_mid[i], vel_mid[i], vel_snd_arr[i]) - DensityFlux(den_mid[i - 1], vel_mid[i - 1], vel_snd_arr[i] ) )
			         +  dt * DensitySource(den_old,vel_old,vel_snd_arr[i]);
			Vel[i] = Vel[i] - (dt / dx) * (VelocityFlux(den_mid[i], vel_mid[i], grad_vel_mid[i], vel_snd_arr[i]) - VelocityFlux(den_mid[i - 1], vel_mid[i - 1], grad_vel_mid[i - 1], vel_snd_arr[i] ) )
					 +  dt * VelocitySource(den_old,vel_old,vel_snd_arr[i]);
			Cur[i] = Vel[i] * Den[i];
		}
} 

/*....................................................................*/	
/*............ Derived Graphene Class  ...............................*/	
/*....................................................................*/	

float GrapheneFluid1D::DensityFlux(float n,float v,float __attribute__((unused)) s){
	float f_1;
	f_1 = n * v;
	return f_1;
}
float GrapheneFluid1D::VelocityFlux(float n,float v,float dv,float s){
	float f_2;
	try {
		if(n<0 || !isfinite(n)){
			throw "ERROR: negative density";
		}
		f_2 = 0.25f * v * v + vel_fer * vel_fer * 0.5f * log(n) + 2.0f * s * s * sqrt(n);//- kin_vis * dv;
	}catch (const char * msg) {
		cerr << msg <<"\nExiting"<<endl;
		exit(EXIT_FAILURE);
	}
	return f_2;
}
float GrapheneFluid1D::DensitySource(__attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0.0f;
}
float GrapheneFluid1D::VelocitySource(__attribute__((unused)) float n,float v,__attribute__((unused)) float s){
	return -1.0f * col_freq * v ;
}


int Fluid1D::GetSnapshotStep() const { return snapshot_step;}
int Fluid1D::GetSnapshotFreq() const {return snapshot_per_period;}


bool Fluid1D::Snapshot() const {
	bool state;
	if(TimeStepCounter % snapshot_step == 0){
		state = true;
	}else{
		state = false;
	}
	return state;
}


void Fluid1D::SaveSnapShot(){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;
	string str_time = to_string(TimeStepCounter / snapshot_step);
	string name_dataset = "snapshot_" + str_time;


	DataSet dataset_den = GrpDen->createDataSet(name_dataset, hdf5_float, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", hdf5_int, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", hdf5_float, atr_dataspace);
	float currenttime=TimeStepCounter * dt;
	atr_step_den.write( hdf5_int, &TimeStepCounter);
	atr_time_den.write( hdf5_float , &currenttime);
	atr_step_den.close();
	atr_time_den.close();
	dataset_den.write(Den, hdf5_float);
	dataset_den.close();

	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, hdf5_float, *DataspaceVelX);
	Attribute atr_step_vel_x = dataset_vel_x.createAttribute("time step", hdf5_int, atr_dataspace);
	Attribute atr_time_vel_x = dataset_vel_x.createAttribute("time", hdf5_float, atr_dataspace);
	atr_step_vel_x.write( hdf5_int, &TimeStepCounter);
	atr_time_vel_x.write( hdf5_float , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();
	dataset_vel_x.write(Vel, hdf5_float);
	dataset_vel_x.close();
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