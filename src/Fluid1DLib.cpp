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
/*
	DenCor = new float[Nx]();
	VelCor = new float[Nx]();
	CurCor = new float[Nx]();
	den_mid = new float[Nx - 1]();
	vel_mid = new float[Nx - 1]();
*/
//	grad_vel_mid = new float[Nx - 1]();
	vel_snd_arr = new float[Nx]();
//	vel_snd_arr_mid = new float[Nx - 1]();
//	lap_den_mid = new float[Nx - 1]();
//	lap_den = new float[Nx]();
///	d3_den_mid = new float[Nx - 1]();
//	d3_den = new float[Nx]();

Umain = new StateVec[Nx]();
Uaux = new StateVec[Nx]();
Umid = new StateVec[Nx-1]();

	SetFDmatrix2(Nx);
	SetFDmatrix3(Nx);
}	

Fluid1D::~Fluid1D() = default;

/*
int Fluid1D::HopscotchFunction(const gsl_vector *x, gsl_vector *f) {
	//auto * params = (struct PhysicalParameters *)p; // ou struct PhysicalParameters * params = (...

	//const float a = params->VSnd;
	//const float b = params->VFer;

//	auto test = static_cast<Fluid1D*>(p);
	const float a = vel_snd;
	const float b = vel_fer;

	//const float a = params->VSnd;
	//const float b = params->VFer;


	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);

	const double y0 = a * (1 - x0);
	const double y1 = b * (x1 - x0 * x0);

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}
*/
/*
int Fluid1D::gslwrapperHopscotchFunction (const gsl_vector *x, void * class_pointer, gsl_vector * f) {
	return static_cast<Fluid1D *>(class_pointer)->HopscotchFunction(x, f);
}
*/

/*
float  Fluid1D::DensityFlux(float n,float v, __attribute__((unused)) float s){
	float f_1;
	f_1 = n * v;
	return f_1;
}
float Fluid1D::VelocityFlux(float n, float v, float dv, float s, float d2n) {
	float f_2;
	f_2 = 0.5f * v * v + n - kin_vis * dv;
	return f_2;
}
*/

/*
float Fluid1D::VelocityFlux(GridPoint1D p, char side) {
	float v= SideAverage(ptr_vel,p,side);
	float n= SideAverage(ptr_den,p,side);
	float dv=SideAverage(ptr_veldx,p,side);
	return 0.5f * v * v + n*vel_snd*vel_snd - kin_vis * dv;
}
*/


float Fluid1D::VelocityFlux(StateVec U) {
	return 0.5f*U.v()*U.v()+U.n()*vel_snd*vel_snd;
}

/*
float Fluid1D::DensityFlux(GridPoint1D p, char side) {
	float v= SideAverage(ptr_vel,p,side);
	float n= SideAverage(ptr_den,p,side);
	return n * v;
}
*/

float Fluid1D::DensityFlux(StateVec U) {
	return U.n()*U.v();
}



float  Fluid1D::DensitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s){
	return 0;
}
float Fluid1D::VelocitySource(float n, float v, float s, float d3den) {
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
		vel_snd_arr[i]= vel_snd;//Sound_Velocity_Anisotropy( static_cast<float>(i)*dx, vel_snd);
	}
//	for(int i = 0; i<Nx-1  ;i++){
//		vel_snd_arr_mid[i]= vel_snd;//Sound_Velocity_Anisotropy( static_cast<float>(i)*dx, vel_snd);
//	}
}
void Fluid1D::SetSound(std::function<float(float)> func) {
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]= func(i*dx);
	}
//	for(int i = 0; i<Nx-1  ;i++){
//		vel_snd_arr_mid[i]= func((i+0.5f)*dx);
//	}
}



void Fluid1D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();

	for (int i = 0; i < Nx; i++ ){
		float noise = (float) rd()/ maxrand ;
	//	Den[i] = 1.0f + 0.005f * (noise - 0.5f);
	//	Vel[i] =0.0f;

		Umain[i].n()= 1.0f + 0.005f * (noise - 0.5f);
		Umain[i].v()=0.0f;
	}
}

void Fluid1D::InitialCondTest(){   //TODO change initial conditions to U stateVec
 	for (int i = 0; i < Nx; i++ ){
		//Vel[i] = 1.0f+tanh(10.0f*(dx*static_cast<float>(i)-0.5f));
		//Den[i]=1.0f+0.05f/(vel_snd*cosh(10.0f*(i*dx-.5f)));
		//Vel[i]=0.0f+0.05f/cosh(10.0f*(i*dx-.5f));
	    //Den[i]=1.0;
		//Vel[i]=  (i<Nx/2) ? 1.0f : -0.5f;
		//Vel[i]=  (i>Nx/3 && i<2*Nx/3 ) ? 1.0f : 0.1f;

		Umain[i].n()=1.0;
	    Umain[i].v()=(i>Nx/3 && i<2*Nx/3 ) ? 1.0f : 0.1f;
	}
}

/*
void Fluid1D::Smooth(int width){
	Average_Filter(Den, DenCor, Nx, width);
	Average_Filter(Vel, VelCor, Nx, width);
	Average_Filter(Cur, CurCor, Nx, width);
}
*/


void Fluid1D::CreateFluidFile(){
	std::string previewfile = "preview_1D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid1D::WriteFluidFile(float t){
	int pos_end = Nx - 1 ;
	int pos_ini = 0;
	if (!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(Vel[pos_end]) ||
	    !isfinite(Vel[pos_ini])) {
		cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
//	data_preview << t << "\t" << Den[pos_end] << "\t" << Vel[pos_end] << "\t" << Den[pos_ini] << "\t" << Vel[pos_ini] << "\n";
	data_preview << t << "\t" << Umain[pos_end/2] << "\n";
}

void Fluid1D::BohmOperator(float bohm) {

	//double beta = 0.001*dt/(dx*dx*dx);
	//gsl_matrix_scale(FDmatrix3, beta); // beta*F3

	gsl_vector *sol = gsl_vector_alloc (Nx);
	gsl_vector *rhs = gsl_vector_alloc (Nx);

	for (int i=0;i<Nx;i++){
		gsl_vector_set(rhs,i,Den[i]); //copiar o array de floats para o vectordouble
	}
	gsl_blas_dgemv(CblasNoTrans, -1.0*bohm*dt/(dx*dx*dx), FDmatrix3, rhs, 0.0, sol);

	for (int i=0;i<Nx;i++){
		Vel[i]=Vel[i]+gsl_vector_get(sol,i);
	}


	//gsl_linalg_LU_solve (BTCSmatrix, permutation_matrix, rhs, sol);
	gsl_vector_free (sol);


}

void Fluid1D::RichtmyerOld(){
	//
	//Calculating the velocity gradient at k time
	//

	//GradientField(Vel, GradVel, dx,  Nx);


	this->RichtmyerStep1();
	//
	//  Calculating the velocity gradient at k+1/2 time
	//
	//GradientField(vel_mid, grad_vel_mid, dx,  Nx-1);

	this->RichtmyerStep2();
//	this->VelocityToCurrent();
}


void Fluid1D::Richtmyer(){
this->RichtmyerStep1();
this->RichtmyerStep2();
}
void Fluid1D::RichtmyerStep1() {
	for ( int i = 0; i <= Nx - 2; i++ ){
		float den_avg   = 0.5f * (Umain[i+1] + Umain[i] ).n();
		float vel_avg   = 0.5f * (Umain[i+1] + Umain[i] ).v();
		Umid[i].n() = den_avg - 0.5f*(dt/dx)*(DensityFlux(Umain[i+1]) - DensityFlux(Umain[i]));
		Umid[i].v() = vel_avg - 0.5f*(dt/dx)*(VelocityFlux(Umain[i+1]) - VelocityFlux(Umain[i]));
	}
}
void Fluid1D::RichtmyerStep2() {
	for ( int i = 1; i <= Nx - 2; i++ ){
		float den_old = Umain[i].n();
		float vel_old = Umain[i].v();
		Umain[i].n() = den_old - (dt/dx)*(DensityFlux(Umid[i]) - DensityFlux(Umid[i-1]));
		Umain[i].v() = vel_old - (dt/dx)*(VelocityFlux(Umid[i]) - VelocityFlux(Umid[i-1]));
	}
}



/*
void Fluid1D::RichtmyerStep1() {
	ChooseGridPointers("MidGrid");
	//
	//Half step calculate density and velocity at time k+0.5 at the spatial midpoints
	//

	for ( int i = 0; i <= Nx - 2; i++ ){
		GridPoint1D midpoint(i, Nx, true);

		den_mid[i] = 0.5f*(Den[i] + Den[i + 1] )
		             - ( 0.5f*dt/dx ) * (DensityFlux(Den[i + 1], Vel[i + 1], vel_snd_arr[i]) - DensityFlux(Den[i], Vel[i], vel_snd_arr[i]) )
		             + ( 0.5f*dt    ) * DensitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
		vel_mid[i] = 0.5f*(Vel[i] + Vel[i + 1] )
		             - ( 0.5f*dt/dx ) * (VelocityFlux(Den[i + 1], Vel[i + 1], GradVel[i + 1], vel_snd_arr[i], lap_den[i+1]) - VelocityFlux(
				Den[i], Vel[i], GradVel[i], vel_snd_arr[i], lap_den[i]))
		             + ( 0.5f*dt    ) * VelocitySource(0.5f * (Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]),
		                                               vel_snd_arr[i], 	0.0 );
	}
}
void Fluid1D::RichtmyerStep2() {
	ChooseGridPointers("MainGrid");
	//
	// Remaining step
	//
	for ( int i = 1; i <= Nx - 2; i++ ){
		GridPoint1D mainpoint(i, Nx, false);
		float den_old = Den[i];
		float vel_old = Vel[i];
		Den[i] = Den[i] - (dt / dx) * (DensityFlux(den_mid[i], vel_mid[i], vel_snd_arr[i]) - DensityFlux(den_mid[i - 1], vel_mid[i - 1], vel_snd_arr[i] ) )
		         +  dt * DensitySource(den_old,vel_old,vel_snd_arr[i]);
		Vel[i] = Vel[i] - (dt / dx) * (VelocityFlux(den_mid[i], vel_mid[i], grad_vel_mid[i], vel_snd_arr[i], lap_den_mid[i]) -
		                               VelocityFlux(
				                               den_mid[i - 1], vel_mid[i - 1], grad_vel_mid[i - 1], vel_snd_arr[i], lap_den_mid[i-1]))
		         +  dt * VelocitySource(den_old, vel_old, vel_snd_arr[i], 0.0f);
	}
}
*/

/*
void Fluid1D::RichtmyerStep1Old() {
	ChooseGridPointers("MidGrid");
	//
	//Half step calculate density and velocity at time k+0.5 at the spatial midpoints
	//

	for ( int i = 0; i <= Nx - 2; i++ ){
		GridPoint1D midpoint(i, Nx, true);

		float den_avg   = 0.5f * (Den[midpoint.E] + Den[midpoint.W] );
		float vel_avg   = 0.5f * (Vel[midpoint.E] + Vel[midpoint.W] );

		den_mid[i] =den_avg - 0.5f*(dt/dx)*(DensityFlux(midpoint, 'E') - DensityFlux(midpoint, 'W'))
		             + ( 0.5f*dt    ) * DensitySource(0.5f*(Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]), vel_snd_arr[i]) ;
		vel_mid[i] = vel_avg - 0.5f*(dt/dx)*(VelocityFlux(midpoint, 'E') - VelocityFlux(midpoint, 'W'))
		             + ( 0.5f*dt    ) * VelocitySource(0.5f * (Den[i] + Den[i + 1]), 0.5f * (Vel[i] + Vel[i + 1]),vel_snd_arr[i], 	0.0 );
	}
}
void Fluid1D::RichtmyerStep2Old() {
	ChooseGridPointers("MainGrid");
	//
	// Remaining step
	//
	for ( int i = 1; i <= Nx - 2; i++ ){
		GridPoint1D mainpoint(i, Nx, false);
		float den_old = Den[i];
		float vel_old = Vel[i];
		Den[i] = den_old - (dt/dx)*(DensityFlux(mainpoint, 'E') - DensityFlux(mainpoint, 'W'))
		         +  dt * DensitySource(den_old,vel_old,vel_snd_arr[i]);
		Vel[i] = vel_old  - (dt/dx)*(VelocityFlux(mainpoint, 'E') - VelocityFlux(mainpoint, 'W'))
		         +  dt * VelocitySource(den_old, vel_old, vel_snd_arr[i], 0.0f);
	}
}*/




void Fluid1D::VelocityToCurrent() {
	for(int i=0; i <= Nx - 1; i++){
		Cur[i] = Vel[i] * Den[i];
	}
}



/*
void Fluid1D::Vliegenthart() {
//WARNING INTENSE NUMERICAL VISCOSITY !!! 
	vector<float>density_new;
	vector<float>velocity_new;

	float den_new, vel_new;

	density_new.push_back(Den[0]);
	velocity_new.push_back(Vel[0]);
	for ( int i = 1; i <= Nx - 2; i++ )
	{
		den_new = 0.5f*(Den[i-1]+Den[i+1])-.5f*(dt / dx) * ( DensityFlux(Den[i+1], Vel[i+1], vel_snd) - DensityFlux(Den[i - 1], Vel[i - 1], vel_snd ) );
		vel_new = 0.5f*(Vel[i-1]+Vel[i+1])-.5f*(dt / dx) * (VelocityFlux(Den[i+1], Vel[i+1],0, vel_snd,0) - VelocityFlux(Den[i - 1], Vel[i - 1],0, vel_snd,0 ) );
	density_new.push_back(den_new);
	velocity_new.push_back(vel_new);
	}
	density_new.push_back(Den[Nx-1]);
	velocity_new.push_back(Vel[Nx-1]);

	for ( int i = 0; i <= Nx - 1; i++ )
	{
		Den[i]=density_new[i];
		Vel[i]=velocity_new[i];
		Cur[i] = Vel[i] * Den[i];
	}
}
*/

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

/*
void Fluid1D::ChooseGridPointers(const string &grid) {
	if(grid == "MidGrid"){  // se ESTÁ na grelha média tem de APONTAR pra outra grelha
		ptr_snd = vel_snd_arr;
		ptr_den = Den;
		ptr_vel = Vel;
		ptr_veldx= GradVel;
		//ptr_dendx = den_dx;
		//ptr_tmp = Tmp;
		//ptr_lap_den = lap_den ;
	}if(grid == "MainGrid"){ // e vice-versa
		ptr_snd = vel_snd_arr_mid;
		ptr_den = den_mid;
		ptr_vel = vel_mid;
		ptr_veldx= grad_vel_mid;

		//ptr_dendx = den_dx_mid;
		//ptr_tmp = tmp_mid;
		//ptr_lap_den = lap_den_mid ;
	}
}
*/
float Fluid1D::SideAverage(const float *input_array, GridPoint1D p, char side) {
	float value;
	switch(side) {
		case 'E': value=input_array[p.E];
			break;
		case 'W': value=input_array[p.W];
			break;
		default: value=0.0f;
	}
	return value;
}

/*
void Fluid1D::SetBTCSmatrix(float beta) {
	BTCSmatrix = gsl_matrix_calloc (Nx, Nx) ;
	gsl_matrix_set_identity(BTCSmatrix);

	gsl_matrix_scale(FDmatrix3, beta); // beta*F3
	gsl_matrix_sub(BTCSmatrix, FDmatrix3); // Identity - (beta*F3)

	permutation_matrix = gsl_permutation_alloc (Nx);
	gsl_linalg_LU_decomp (BTCSmatrix, permutation_matrix, &permutation_index_s);

}
*/

void Fluid1D::RungeKuttaTVD() {
	float DenNumFluxW;
	float DenNumFluxE;
	float VelNumFluxW;
	float VelNumFluxE;
	StateVec UEleft(Umain[0]);
	StateVec UEright(Umain[0]);
	StateVec UWleft(Umain[0]);
	StateVec UWright(Umain[0]);

	for (int i = 1; i < Nx-1; ++i) {

		CellHandler1D cell(i, this, Umain);

		UEleft  = cell.TVD('E','L');
		UEright = cell.TVD('E','R');
		UWleft  = cell.TVD('W','L');
		UWright = cell.TVD('W','R');

		DenNumFluxW= NumericalFlux::Central(this,UWleft,UWright).n();
		DenNumFluxE= NumericalFlux::Central(this,UEleft,UEright).n();
		VelNumFluxW= NumericalFlux::Central(this,UWleft,UWright).v();
		VelNumFluxE= NumericalFlux::Central(this,UEleft,UEright).v();

		Uaux[i].n()=Umain[i].n()-(dt/dx)*(DenNumFluxW-DenNumFluxE);
		Uaux[i].v()=Umain[i].v()-(dt/dx)*(VelNumFluxW-VelNumFluxE);
	}
	for (int i = 1; i < Nx-1; ++i) {
		CellHandler1D cell(i, this, Uaux);
		UEleft  = cell.TVD('E','L');
		UEright = cell.TVD('E','R');
		UWleft  = cell.TVD('W','L');
		UWright = cell.TVD('W','R');

		DenNumFluxW= NumericalFlux::Central(this,UWleft,UWright).n();
		DenNumFluxE= NumericalFlux::Central(this,UEleft,UEright).n();
		VelNumFluxW= NumericalFlux::Central(this,UWleft,UWright).v();
		VelNumFluxE= NumericalFlux::Central(this,UEleft,UEright).v();

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

float Fluid1D::JacobianSignum(StateVec U, std::string key) {

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
