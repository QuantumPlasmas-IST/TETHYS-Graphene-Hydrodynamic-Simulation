#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>




#include "TethysLib.h"
#include "Tethys1DLib.h"
#include <H5Cpp.h>

using namespace H5;
using namespace std;


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905
#endif

#ifndef C_SPEED
#	define C_SPEED 1000.0
#endif

GrapheneFluid1D::GrapheneFluid1D(int size_n, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency): Fluid1D(size_n, sound_velocity, shear_viscosity){
	vel_fer =fermi_velocity;
	col_freq =collision_frequency;
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}

/*....................................................................*/	
/*.......... 1 Dimensional Fluid Class ...............................*/	
/*....................................................................*/	
Fluid1D::Fluid1D(int size_nx, float sound_velocity, float shear_viscosity): TETHYSBase{size_nx, 0, 1}{
	Nx = size_nx;
	vel_snd =sound_velocity;
	kin_vis =shear_viscosity;
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
float  Fluid1D::DensitySource(float n,float v,float s){
	return 0;
}
float  Fluid1D::VelocitySource(float n,float v,float s){
	return 0;
}

void Fluid1D::CFLCondition(){
		dx = leng / ( float ) ( Nx - 1 );
		dt = dx/10.0f;
}

void Fluid1D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx-1  ;i++){
		vel_snd_arr[i]= Sound_Velocity_Anisotropy(i, dx, vel_snd);
	}
}
		
void Fluid1D::InitialCondRand(){
	srand (static_cast<unsigned int>(time(NULL)));
	for (int i = 0; i < Nx; i++ ){
		float noise = (float) rand()/ (float) RAND_MAX ;
		Den[i] = 1.0f + 0.005f * (noise - 0.5f);
	}
}

void Fluid1D::InitialCondTest(){
	float mean=dx*Nx/2.0f;
	float sigma=dx*Nx/10.0f;
	for (int i = 0; i < Nx; i++ ){
		//Vel[i] = 0.4f*exp(-0.5f*((i*dx-mean)*(i*dx-mean)/(sigma*sigma)))/sigma;
		Vel[i] = 1.0f+tanh(10.0f*(dx*i-0.5));
	}
}

void Fluid1D::SetKinVis(float x){ kin_vis=x;}
void Fluid1D::SetVelSnd(float x){ vel_snd=x; }
float Fluid1D::GetVelSnd(){ return vel_snd; }
float Fluid1D::GetKinVis(){ return kin_vis; }
float Fluid1D::GetDx(){return dx;}
float Fluid1D::GetDt(){return dt;}


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
	int pos_end = Nx - 1 ;
	int pos_ini = 0;
	try {
		if (!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(Vel[pos_end]) ||
		    !isfinite(Vel[pos_ini])) {
			throw "ERROR: numerical method failed to converge";
		}
//data_preview << t << "\t" << DenCor[Nx - 1] << "\t" << VelCor[Nx - 1] << "\t" << DenCor[0] << "\t" << VelCor[0] << "\n";
		data_preview << t << "\t" << Den[pos_end] << "\t" << Vel[pos_end] << "\t" << Den[pos_ini] << "\t" << Vel[pos_ini] << "\n";
	}catch (const char* msg) {
		cerr << msg << endl;
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
					;//+ ( 0.5f*dt    ) * DensitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
			vel_mid[i] = 0.5f*(Vel[i] + Vel[i + 1] )
				- ( 0.5f*dt/dx ) * (VelocityFlux(Den[i + 1], Vel[i + 1], GradVel[i + 1], vel_snd_arr[i]) - VelocityFlux(Den[i], Vel[i], GradVel[i], vel_snd) )
					;//+ ( 0.5f*dt    ) * VelocitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
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
			Den[i] = Den[i] - (dt / dx) * (DensityFlux(den_mid[i], vel_mid[i], vel_snd_arr[i]) - DensityFlux(den_mid[i - 1], vel_mid[i - 1], vel_snd) )
			        ;// +  dt * DensitySource(den_old,vel_old,vel_snd_arr[i]);
			Vel[i] = Vel[i] - (dt / dx) * (VelocityFlux(den_mid[i], vel_mid[i], grad_vel_mid[i], vel_snd_arr[i]) - VelocityFlux(den_mid[i - 1], vel_mid[i - 1], grad_vel_mid[i - 1], vel_snd) )
			        ;// +  dt * VelocitySource(den_old,vel_old,vel_snd_arr[i]);
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
	} catch (const char * msg) {
		cerr << msg <<endl;
		exit(EXIT_FAILURE);
	}
	return f_2;
}
float GrapheneFluid1D::DensitySource(float n,float v,float s){
	float q_1=0.0f;
	return q_1;
}
float GrapheneFluid1D::VelocitySource(float n,float v,float s){
	float q_2;
	q_2= -1.0f * col_freq * (v - 1);
	return q_2;
}

void GrapheneFluid1D::SetVelFer(float x){ vel_fer=x; }
float GrapheneFluid1D::GetVelFer(){ return vel_fer; }
void GrapheneFluid1D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid1D::GetColFreq(){ return col_freq; }

void Fluid1D::SaveSnapShot(int time_step,int snapshot_step){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);

	string str_time = to_string(time_step / snapshot_step);
	string name_dataset = "snapshot_" + str_time;

	DataSet dataset_den = GrpDen->createDataSet(name_dataset, hdf5_float, *DataspaceDen);
	dataset_den.write(Den, hdf5_float);
	dataset_den.close();

	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, hdf5_float, *DataspaceVelX);
	dataset_vel_x.write(Vel, hdf5_float);
	dataset_vel_x.close();

}


void GrapheneFluid1D::WriteAtributes(){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	int total_steps= static_cast<int>(Tmax / dt);
	//Create the data space for the attribute.
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = GrpDat->createAttribute("S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = GrpDat->createAttribute("Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_kin_vis = GrpDat->createAttribute("Kinetic viscosity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = GrpDat->createAttribute("Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = GrpDat->createAttribute("Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = GrpDat->createAttribute("Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = GrpDat->createAttribute("Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_space_points = GrpDat->createAttribute("Number of spatial points", hdf5_int, atr_dataspace);
	Attribute atr_num_time_steps = GrpDat->createAttribute("Number of time steps", hdf5_int, atr_dataspace);
	// Write the attribute data.
	atr_vel_snd.write( hdf5_float, &vel_snd);
	atr_vel_fer.write( hdf5_float, &vel_fer);
	atr_col_freq.write(hdf5_float, &col_freq);
	atr_kin_vis.write(hdf5_float, &kin_vis); 
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_num_space_points.write( hdf5_int, &Nx);
	atr_total_time.write( hdf5_float, &Tmax);
	atr_num_time_steps.write(hdf5_int, &total_steps);
	// Close the attributes.
	atr_num_time_steps.close();
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
}




void GrapheneFluid1D::CFLCondition(){
	dx = leng / ( float ) ( Nx - 1 );	
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
}	



/*....................................................................*/

void ElectroAnalysis::CreateElectroFile(GrapheneFluid1D& graphene){
	//graphene.SetFileName();
	std::string infix = graphene.GetInfix();
	std::string electrofile = "electro_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific; 	
}

void ElectroAnalysis::WriteElectroFile(float t,GrapheneFluid1D& graphene){
		float q_net = this->NetCharge(graphene);
		float i_avg = this->AverageCurrent(graphene);
		float p_ohm = this->OhmPower(graphene);
		float dipole_var=this->ElectricDipoleVariation(graphene);
		float dipole=this->ElectricDipole(graphene);
		data_electro << t << "\t" << q_net << "\t" << i_avg << "\t" << q_net * q_net * 0.5 << "\t" << p_ohm << "\t" << dipole << "\t" << dipole_var << "\n";
}

float ElectroAnalysis::NetCharge(GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.DenCor);
}

float ElectroAnalysis::AverageCurrent(GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.CurCor);
}

float ElectroAnalysis::ElectricDipoleVariation(GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.CurCor);
}

float ElectroAnalysis::ElectricDipole(GrapheneFluid1D& graphene){
	float p=0.0;
	float dx=graphene.GetDx();
	for(int j=1;j<graphene.SizeX()/2;j++){	
		p += dx*(2*j-2)*graphene.DenCor[2 * j - 2] + 4 * dx * (2 * j - 1) * graphene.DenCor[2 * j - 1] + dx * (2 * j) * graphene.DenCor[2 * j];
	}
	p = p*graphene.GetDx()/3.0f;
	return p;
}

float ElectroAnalysis::OhmPower(GrapheneFluid1D& graphene){
	float itg=0.0;
	for(int j=1;j<graphene.SizeX()/2;j++){
		itg += graphene.CurCor[2 * j - 2] * graphene.VelCor[2 * j - 2] + 4 * graphene.CurCor[2 * j - 1] * graphene.VelCor[2 * j - 1] + graphene.CurCor[2 * j] * graphene.VelCor[2 * j];
	}
	itg = itg*graphene.GetDx()/3.0f;
	return itg;	
}

