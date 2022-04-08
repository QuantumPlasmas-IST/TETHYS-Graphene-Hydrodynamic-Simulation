/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/Fluid1DLib.h"
#include "includes/Cell1DLib.h"
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
	GradVel= new float[Nx]();
	Cur = new float[Nx]();
	vel_snd_arr = new float[Nx]();

	Umain = new StateVec[Nx]();
	Uaux = new StateVec[Nx]();
	Umid = new StateVec[Nx-1]();

}	

Fluid1D::~Fluid1D() = default;


float Fluid1D::VelocityFlux( StateVec U) {
	return 0.5f*U.v()*U.v()+U.n()*U.S()*U.S();
}



float Fluid1D::DensityFlux(StateVec U) {
	return U.n()*U.v();
}

float  Fluid1D::DensitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0;
}
float Fluid1D::VelocitySource(float n, float v, float s, float d3den) {
	return 0;
}

float  Fluid1D::DensitySource(StateVec U){
	return 0;
}
float Fluid1D::VelocitySource(StateVec U) {
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
		Umain[i].n()=1.0;
	    Umain[i].v()=(i>Nx/3 && i<2*Nx/3 ) ? 1.0f : 0.1f;
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
	RichtmyerStep1();
	RichtmyerStep2();
}
void Fluid1D::RichtmyerStep1() {
	for ( int i = 0; i <= Nx - 2; i++ ){
		float den_avg   = 0.5f * (Umain[i+1] + Umain[i] ).n();
		float vel_avg   = 0.5f * (Umain[i+1] + Umain[i] ).v();
		Umid[i].n() = den_avg - 0.5f*(dt/dx)*(DensityFlux(Umain[i+1]) - DensityFlux(Umain[i]))
							+ (0.5f*dt) * DensitySource(0.5f * (Umain[i+1] + Umain[i] )) ;
		Umid[i].v() = vel_avg - 0.5f*(dt/dx)*(VelocityFlux(Umain[i+1]) - VelocityFlux(Umain[i]))
							+ (0.5f*dt) * VelocitySource(0.5f * (Umain[i+1] + Umain[i] )) ;
	}
}
void Fluid1D::RichtmyerStep2() {
	for ( int i = 1; i <= Nx - 2; i++ ){
		StateVec Uold(Umain[i]);
		float den_old = Uold.n();
		float vel_old = Uold.v();
//		float den_old = Umain[i].n();
//		float vel_old = Umain[i].v();
		Umain[i].n() = den_old - (dt/dx)*(DensityFlux(Umid[i]) - DensityFlux(Umid[i-1]))
				                 + dt * DensitySource(Uold);
		Umain[i].v() = vel_old - (dt/dx)*(VelocityFlux(Umid[i]) - VelocityFlux(Umid[i-1]))
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
	snapshot_step = 1; //points_per_period / snapshot_per_period;

	string str_time = to_string(TimeStepCounter );/// snapshot_step);
	str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
	string name_dataset = "snapshot_" + str_time;

	DataSet dataset_den = GrpDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceDen);
	Attribute atr_step_den = dataset_den.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_den = dataset_den.createAttribute("time", HDF5FLOAT, atr_dataspace);
	float currenttime= static_cast<float>(TimeStepCounter) * dt;
	atr_step_den.write(HDF5INT, &TimeStepCounter);
	atr_time_den.write(HDF5FLOAT , &currenttime);
	atr_step_den.close();
	atr_time_den.close();
	dataset_den.write(Den, HDF5FLOAT);
	dataset_den.close();

	DataSet dataset_vel_x = GrpVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceVelX);
	Attribute atr_step_vel_x = dataset_vel_x.createAttribute("time step", HDF5INT, atr_dataspace);
	Attribute atr_time_vel_x = dataset_vel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
	atr_step_vel_x.write(HDF5INT, &TimeStepCounter);
	atr_time_vel_x.write(HDF5FLOAT , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();
	dataset_vel_x.write(Vel, HDF5FLOAT);
	dataset_vel_x.close();
}

void Fluid1D::RungeKuttaTVD() {
	float DenNumFluxW;
	float DenNumFluxE;
	float VelNumFluxW;
	float VelNumFluxE;
	StateVec UEleft{};
	StateVec UEright{};
	StateVec UWleft{};
	StateVec UWright{};

	for (int i = 1; i < Nx-1; ++i) {

		CellHandler1D cell(i, this, Umain);
	//	UEleft  = cell.TVD(Umain,i,'E','L');
	//	UEright = cell.TVD(Umain,i,'E','R');
	//	UWleft  = cell.TVD(Umain,i,'W','L');
	//	UWright = cell.TVD(Umain,i,'W','R');

		// redefini os StateVec's, para estarem iguais à literatura para eu me orientar melhor :(
		UEleft  = Umain[i];
		UEright = Umain[i-1];
		UWleft  = Umain[i+1];
		UWright = Umain[i];

		// calculei os laplacianos da raiz quadrada da densidade usando o método das diferenças centradas
		float lap_sqrtn = (Umain[i+1].n()-2*Umain[i].n()+Umain[i-1].n())/(dx*dx);
		UEleft.d2den()  = lap_sqrtn;
		UEright.d2den() = lap_sqrtn;
		UWleft.d2den()  = lap_sqrtn;
		UWright.d2den() = lap_sqrtn;

		// alterei o método para a média para começar
		DenNumFluxE = NumericalFlux::Central(this,UEleft,UEright).n();
		DenNumFluxW = NumericalFlux::Central(this,UWleft,UWright).n();
		VelNumFluxE = NumericalFlux::Central(this,UEleft,UEright).v();
		VelNumFluxW = NumericalFlux::Central(this,UWleft,UWright).v();

		Uaux[i].n()=Umain[i].n()-(dt/dx)*(DenNumFluxW-DenNumFluxE);
		Uaux[i].v()=Umain[i].v()-(dt/dx)*(VelNumFluxW-VelNumFluxE);
	}
	for (int i = 1; i < Nx-1; ++i) {
		CellHandler1D cell(i, this, Uaux); // vetor Uaux vazio?

	//	UEleft  = cell.TVD(Uaux,i,'E','L');
	//	UEright = cell.TVD(Uaux,i,'E','R');
	//	UWleft  = cell.TVD(Uaux,i,'W','L');
	//	UWright = cell.TVD(Uaux,i,'W','R');

		// redefini os StateVec's, para estarem iguais à literatura para eu me orientar melhor :(
		UEleft  = Uaux[i];
		UEright = Uaux[i-1];
		UWleft  = Uaux[i+1];
		UWright = Uaux[i];

		// calculei os laplacianos da raiz quadrada da densidade usando o método das diferenças centradas
		float lap_sqrtn = (Umain[i+1].n()-2*Umain[i].n()+Umain[i-1].n())/(dx*dx);
		UEleft.d2den()  = lap_sqrtn;
		UEright.d2den() = lap_sqrtn;
		UWleft.d2den()  = lap_sqrtn;
		UWright.d2den() = lap_sqrtn;

		DenNumFluxE= NumericalFlux::Central(this,UEleft,UEright).n();
		DenNumFluxW= NumericalFlux::Central(this,UWleft,UWright).n();
		VelNumFluxE= NumericalFlux::Central(this,UEleft,UEright).v();
		VelNumFluxW= NumericalFlux::Central(this,UWleft,UWright).v();

		Umain[i].n()=0.5f*(Umain[i].n()+Uaux[i].n())-(0.5f*dt/dx)*(DenNumFluxW-DenNumFluxE);
		Umain[i].v()=0.5f*(Umain[i].v()+Uaux[i].v())-(0.5f*dt/dx)*(VelNumFluxW-VelNumFluxE);
	}
}

void Fluid1D::McCormack() {
	for (int i = 1; i < Nx-1; ++i) {
		Uaux[i].n()=Umain[i].n()-(dt/dx)*(DensityFlux(Umain[i+1])-DensityFlux(Umain[i]));
		Uaux[i].v()=Umain[i].v()-(dt/dx)*(VelocityFlux(Umain[i+1])-VelocityFlux(Umain[i]));
	}
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=0.5f*(Umain[i].n()+Uaux[i].n())-(0.5f*dt/dx)*(DensityFlux(Uaux[i])-DensityFlux(Uaux[i-1]));
		Umain[i].v()=0.5f*(Umain[i].v()+Uaux[i].v())-(0.5f*dt/dx)*(VelocityFlux(Uaux[i])-VelocityFlux(Uaux[i-1]));
	}
}

void Fluid1D::Upwind(){
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=Umain[i].n()-(dt/dx)*(DensityFlux(Umain[i])-DensityFlux(Umain[i-1]));
		Umain[i].v()=Umain[i].v()-(dt/dx)*(VelocityFlux(Umain[i])-VelocityFlux(Umain[i-1]));
	}
}

void Fluid1D::LaxFriedrichs(){
	for (int i = 1; i < Nx-1; ++i) {
		Umain[i].n()=Den[i];
		Umain[i].v()=Vel[i];
		Den[i]=0.5f*(Umain[i-1].n()+Umain[i+1].n())-0.5f*(dt/dx)*(DensityFlux(Umain[i+1])-DensityFlux(Umain[i-1]));
		Vel[i]=0.5f*(Umain[i-1].v()+Umain[i+1].v())-0.5f*(dt/dx)*(VelocityFlux(Umain[i+1])-VelocityFlux(Umain[i-1]));
	}
}

float Fluid1D::JacobianSpectralRadius(StateVec U) {
	float l1=abs(U.v()+ vel_snd*sqrt(U.n()));
	float l2=abs(U.v()- vel_snd*sqrt(U.n()));
	return max(l1,l2);
}

StateVec Fluid1D::ConservedFlux(StateVec U) {
	StateVec Uout{};
	Uout.n()= this->DensityFlux(U);
	Uout.v()= this->VelocityFlux(U);
	return Uout;
}

float Fluid1D::JacobianSignum( StateVec U, std::string key) {

	float l1= Signum(U.v()+vel_snd*sqrt(U.n()));
	float l2= Signum(U.v()-vel_snd*sqrt(U.n()));
	float entry=0;
	if(key=="11"){
		entry=l1*0.5f;
	}else if(key=="12"){
		entry=l2* sqrt(U.n())/(2.0f*vel_snd);
	}else if(key=="21"){
		entry=l1* vel_snd/sqrt(U.n());
	}else if(key=="22"){
		entry=l2*0.5f;
	}else entry=0.0f;
return entry;
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

