/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/Fluid1DLib.h"
//#include "includes/Cell1DLibNOT.h"
#include "includes/SetUpParametersLib.h"



using namespace H5;
using namespace std;


Fluid1D::Fluid1D(const SetUpParameters &input_parameters) : TethysBase{input_parameters.SizeX, 0, 1}{
	Nx = input_parameters.SizeX;
	vel_snd = input_parameters.SoundVelocity;
	kin_vis = input_parameters.ShearViscosity;
	col_freq=input_parameters.CollisionFrequency;

	param = {vel_snd,0.0f,0.0f,kin_vis,0.0f,0.0f,col_freq,0.0f};

	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	Den = new float[Nx]();
	Vel = new float[Nx]();
	//GradVel= new float[Nx]();
	Cur = new float[Nx]();
	vel_snd_arr = new float[Nx]();

	Umain = new StateVec1D[Nx]();
	Uaux = new StateVec1D[Nx]();
	Umid = new StateVec1D[Nx - 1]();

}	

Fluid1D::~Fluid1D() = default;


float Fluid1D::VelocityFlux(StateVec1D U) {
	return 0.5f*U.v()*U.v() + vel_fer * vel_fer  *0.5f* log(U.n()+1.0E-6f) ;//- kin_vis*U.grad_v();
}



float Fluid1D::DensityFlux(StateVec1D U) {
	return U.n()*U.v();
}

float  Fluid1D::DensitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0;
}
float Fluid1D::VelocitySource(float n, float v, float s, float d3den) {
	return 0;
}

float  Fluid1D::DensitySource(StateVec1D U){
	return 0;
}
float Fluid1D::VelocitySource(StateVec1D U) {
	return  -1.0f * col_freq * U.v() ;
}


void Fluid1D::CflCondition(){
		dx = lengX / ( float ) ( Nx - 1 );
		dt = dx/10.0f;
}

void Fluid1D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]= vel_snd;
		Umain[i].S()=vel_snd;
		Uaux[i].S()=vel_snd;
	}
	for(int i = 0; i<Nx-1  ;i++){
		Umid[i].S()=vel_snd;
	}
}
void Fluid1D::SetSound(const std::function<float(float)>& func) {
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]= func(i*dx);
		Umain[i].S()=func(i*dx);
		Uaux[i].S()=func(i*dx);
	}
	for(int i = 0; i<Nx-1  ;i++){
		Umid[i].S()=func((i+0.5f)*dx);
	}
}



void Fluid1D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();

	for (int i = 0; i < Nx; i++ ){
		float noise = (float) rd()/ maxrand ;
		Umain[i].n()= 1.0f + 0.0001f * (noise - 0.5f);
		Umain[i].v()= 0.0f;
	}
	this->SetSound();
}

void Fluid1D::InitialCondTest(){
 	for (int i = 0; i < Nx; i++ ){
		//Umain[i].v()= 1.5f; // (i>3*Nx/8 && i<5*Nx/8 ) ? 3.0f : 0.0f; //1.5f;//
	    Umain[i].v()= 1.0f/(1.0f+5.0f* pow(cosh((i*dx-0.5f)*12.0f),2.f));
	    Umain[i].n()= 0.2f+0.2f/ pow(cosh((i*dx-0.5f)*12.0f),2.f); //(i>3*Nx/8 && i<5*Nx/8 ) ? 1.0f : 0.1f; //0.2f+0.2f/ pow(cosh((i*dx-0.5f)*12.0f),2);//
	}
	this->SetSound();
}
void Fluid1D::InitialCondGeneral(function<float(float)> fden, function<float(float)> fvx) {
	float x;
	for (int i = 0; i < Nx; ++i) {
		x=i*dx;
		Umain[i].n()=fden(x);
		Umain[i].v()=fvx(x);
	}
	this->SetSound();
}



void Fluid1D::CreateFluidFile(){
	std::string previewfile = "preview_1D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid1D::WriteFluidFile(float t){
	int pos_end = Nx - 1 ;
	int pos_ini = 0;
	if (!isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].v()) ||
	    !isfinite(Umain[pos_ini].v())) {
		cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
//data_preview << t << "\t" << Umain[pos_end/2] << "\n";
	data_preview << t << "\t" << Umain[pos_ini] << "\t" << Umain[pos_end] << "\n";
}


void Fluid1D::Richtmyer(){
	if(kin_vis!=0){
		CalcVelocityGradient(Umain,Nx);
	}
	RichtmyerStep1();
	if(kin_vis!=0) {
		CalcVelocityGradient(Umid, Nx - 1);
	}
	RichtmyerStep2();
}
void Fluid1D::RichtmyerStep1() {
	for ( int i = 0; i <= Nx - 2; i++ ){
//		float den_avg   = 0.5f * (Umain[i+1] + Umain[i] ).n();
//		float vel_avg   = 0.5f * (Umain[i+1] + Umain[i] ).v();

		StateVec1D Uavg{};
		Uavg = 0.5f*(Umain[i+1] + Umain[i]);

		Umid[i].n() = Uavg.n() - 0.5f*(dt/dx)*(DensityFlux(Umain[i+1]) - DensityFlux(Umain[i]))
							+ (0.5f*dt) * DensitySource(Uavg) ;
		Umid[i].v() = Uavg.v() - 0.5f*(dt/dx)*(VelocityFlux(Umain[i+1]) - VelocityFlux(Umain[i]))
							+ (0.5f*dt) * VelocitySource(Uavg) ;
	}
}
void Fluid1D::RichtmyerStep2() {
	for ( int i = 1; i <= Nx - 2; i++ ){
		StateVec1D Uold(Umain[i]);
//		float den_old = Uold.n();
//		float vel_old = Uold.v();
		Umain[i].n() = Uold.n() - (dt/dx)*(DensityFlux(Umid[i]) - DensityFlux(Umid[i-1]))
				                 + dt * DensitySource(Uold);
		Umain[i].v() = Uold.v() - (dt/dx)*(VelocityFlux(Umid[i]) - VelocityFlux(Umid[i-1]))
		                         + dt * VelocitySource(Uold);
	}
}




void Fluid1D::VelocityToCurrent() {
	for(int i=0; i <= Nx - 1; i++){
		Cur[i] = Vel[i] * Den[i];
	}
}

int Fluid1D::GetSnapshotStep() const { return snapshot_step;}
int Fluid1D::GetSnapshotFreq() const {return snapshot_per_period;}


bool Fluid1D::Snapshot() const {
	return TimeStepCounter % snapshot_step == 0;
}


void Fluid1D::SaveSnapShot(){

	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;

    this->CopyFields();

	string str_time = to_string(TimeStepCounter / snapshot_step );/// snapshot_step);
	//TODO TRATAR AQUI DISTO ta a dar asneira porque o numero dos snapshotsé muito pequeno parece quando chega ao snapshot 100000
	str_time.insert(str_time.begin(), 7 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

    float currenttime= static_cast<float>(TimeStepCounter) * dt;

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
    dataset_vel_x.write(Vel, HDF5FLOAT);
    dataset_vel_x.close();
    atr_step_vel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_x.write(HDF5FLOAT , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();


/*
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
    */

}


StateVec1D Fluid1D::ConservedFlux(StateVec1D U) {
	StateVec1D Uout{};
	Uout.n()= this->DensityFlux(U);
	Uout.v()= this->VelocityFlux(U);
	return Uout;
}


void Fluid1D::CopyFields() {
	for (int i = 0; i < Nx; ++i) {
		Den[i]=Umain[i].n();
		Vel[i]=Umain[i].v();
		Cur[i]=Umain[i].v()*Umain[i].n();
	}
}


void Fluid1D::SaveSound() {
	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(vel_snd_arr, HDF5FLOAT);
	dataset_vel_snd.close();
}

void Fluid1D::CalcVelocityGradient(StateVec1D * u_vec, int size_x) {
	for ( int i = 1; i < size_x-1 ; i++ )
	{
		u_vec[i].grad_v() = (-0.5f * u_vec[i - 1].v() + 0.5f * u_vec[i + 1].v()) / dx;
	}
	u_vec[0].grad_v() = (-1.5f * u_vec[0].v() + 2.0f * u_vec[1].v() - 0.5f * u_vec[2].v()) / dx;
	u_vec[size_x - 1].grad_v() = (0.5f * u_vec[size_x - 1 - 2].v() - 2.0f * u_vec[size_x - 1 - 1].v() + 1.5f * u_vec[size_x - 1].v()) / dx;
}

void Fluid1D::CalcVelocityLaplacian(StateVec1D * u_vec, int size_x) {
	for ( int i = 1; i < size_x-1 ; i++ )
	{
		u_vec[i].grad_v() = (u_vec[i - 1].v() -2.0f*u_vec[i].v() + u_vec[i + 1].v()) / (dx*dx);
	}
	u_vec[0].grad_v() = (2.0f * u_vec[0].v() - 5.0f * u_vec[1].v() + 4.0f * u_vec[2].v()-u_vec[3].v()) / (dx*dx);
	u_vec[size_x - 1].grad_v() = ( 2.0f * u_vec[size_x - 1].v() - 5.0f * u_vec[size_x - 1-1].v() + 4.0f * u_vec[size_x - 1-2].v()-u_vec[size_x - 1-3].v()) / (dx*dx);
}

void Fluid1D::ParabolicFTCS() {
	CalcVelocityLaplacian(Umain,Nx);

	for (int i = 0; i < Nx ; ++i) {
		Umain[i].v() = Umain[i].v() + kin_vis * dt *Umain[i].grad_v() ; //TODO mudar o nome para lap_v se isto funcionar
	}
}
