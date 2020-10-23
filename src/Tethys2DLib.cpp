// 2D version



#include "Tethys2DLib.h"


using namespace H5;
using namespace std;


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif

Fluid2D::Fluid2D(int size_nx, int size_ny, const SetUpInput &input_parameters) : TethysBase{size_nx, size_ny, 2}{
	Nx = size_nx;
	Ny = size_ny;

	vel_snd = input_parameters.SoundVelocity;//sound_velocity;
	kin_vis = input_parameters.ShearViscosity;//shear_viscosity;

	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	// main grid variables Nx*Ny
	Den 		= new float[Nx * Ny]();
	VelX 		= new float[Nx * Ny]();
	VelY 		= new float[Nx * Ny]();
	FlxX 		= new float[Nx * Ny]();
	FlxY 		= new float[Nx * Ny]();
	CurX 		= new float[Nx * Ny]();
	CurY 		= new float[Nx * Ny]();
	vel_snd_arr	= new float[Nx*Ny]();

	lap_flxX = new float[Nx*Ny](); //new grids for the laplacians
	lap_flxY = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();  
	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();
	vel_snd_arr_mid	= new float[(Nx-1)*(Ny-1)]();
}
	
Fluid2D::~Fluid2D(){
	delete [] Den;
	delete [] VelX;
	delete [] VelY;
	delete [] FlxX;
	delete [] FlxY;
	delete [] CurX;
	delete [] CurY;
	delete [] den_mid;
	delete [] flxX_mid;
	delete [] flxY_mid;
	delete [] lap_flxX;
	delete [] lap_flxY;
	delete [] vel_snd_arr;
}



void Fluid2D::SetSound(){
	for(int kp=0; kp<=Nx*Ny-1; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		int j = divresult.quot;
		int i = divresult.rem;
		vel_snd_arr[kp]= Sound_Velocity_Anisotropy(i*dx, j*dy , vel_snd);
	}
	for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++) { //correr todos os pontos da grelha secundaria
		div_t divresult;
		divresult = div(ks, Nx - 1);
		int j = divresult.quot;
		int i = divresult.rem;
		vel_snd_arr_mid[ks]= Sound_Velocity_Anisotropy((i+0.5f)*dx, (j+0.5f)*dy , vel_snd);
	}


}




void Fluid2D::InitialCondRand(){
	random_device rd;
	float maxrand;
	maxrand = (float) rd.max();
	//srand (static_cast<unsigned int>(time(NULL)));

	for (int i = 0; i < Nx; i++ ){
		for (int j=0; j<Ny; j++){
		float noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
			Den[i + j * Nx] = 1.0f + 0.005f * (noise - 0.5f);
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


void Fluid2D::MassFluxToVelocity(){
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] / Den[c];
		VelY[c]= FlxY[c] / Den[c];
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}

void Fluid2D::Richtmyer(){
		int northeast,northwest,southeast,southwest;
		float den_north, den_south ,den_east ,den_west, px_north, px_south, px_east, px_west, py_north, py_south, py_east, py_west,m_east,m_west,m_north,m_south;
		float  sound_north, sound_south ,sound_east ,sound_west;
		//k=i+j*Nx
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
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

			//posso definir aqui a "massa" nos 4 ponto s cardeais
			m_east=pow(den_east,1.5f); // e assim sucessivamente m_west m_north m_south que depois sao reutilizadeas nos 12 fluxos
			m_west=pow(den_west,1.5f);
			m_north=pow(den_north,1.5f);
			m_south=pow(den_south,1.5f);

			float DenAvg = 0.25f*(Den[southwest] + Den[southeast] + Den[northwest] + Den[northeast]);
			float FlxXAvg = 0.25f*(FlxX[southwest] + FlxX[southeast] + FlxX[northwest] + FlxX[northeast]);
			float FlxYAvg = 0.25f*(FlxY[southwest] + FlxY[southeast] + FlxY[northwest] + FlxY[northeast]);
			den_mid[ks] = DenAvg
							-0.5f*(dt/dx)*(
								DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
								DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
							-0.5f*(dt/dy)*(
								DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
								DensityFluxY(den_south, px_south, py_south,m_south,sound_south))//;
							+0.5f*dt*DensitySource(DenAvg,FlxXAvg,FlxYAvg,0.0f,0.0f);
			flxX_mid[ks] = FlxXAvg
							-0.5f*(dt/dx)*(
								MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
								MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
							-0.5f*(dt/dy)*(
								MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
								MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))//;
							+0.5f*dt*MassFluxXSource(DenAvg,FlxXAvg,FlxYAvg,0.0f,0.0f);
			flxY_mid[ks] = FlxYAvg
							-0.5f*(dt/dx)*(
								MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
								MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
							-0.5f*(dt/dy)*(
								MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
								MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))//;
							+0.5f*dt*MassFluxYSource(DenAvg,FlxXAvg,FlxYAvg,0.0f,0.0f);
		}
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
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

				//posso definir aqui a "massa" nos 4 ponto s cardeais
				m_east=pow(den_east,1.5f); // e assim sucessivamente m_west m_north m_south que depois sao reutilizadeas nos 12 fluxos
				m_west=pow(den_west,1.5f);
				m_north=pow(den_north,1.5f);
				m_south=pow(den_south,1.5f);
				float DenOld = Den[kp];
				float FlxXOld = FlxX[kp];
				float FlxYOld = FlxY[kp];
				Den[kp] = DenOld-(dt/dx)*(
									DensityFluxX(den_east, px_east, py_east,m_east,sound_east)-
									DensityFluxX(den_west, px_west, py_west,m_west,sound_west))
									-(dt/dy)*(
									DensityFluxY(den_north, px_north, py_north,m_north,sound_north)-
									DensityFluxY(den_south, px_south, py_south,m_south,sound_south))//;
									+dt*DensitySource(DenOld,FlxXOld,FlxYOld,0.0f,0.0f);
				FlxX[kp] = FlxXOld-(dt/dx)*(
									MassFluxXFluxX(den_east, px_east, py_east,m_east,sound_east)-
									MassFluxXFluxX(den_west, px_west, py_west,m_west,sound_west))
									-(dt/dy)*(
									MassFluxXFluxY(den_north, px_north, py_north,m_north,sound_north)-
									MassFluxXFluxY(den_south, px_south, py_south,m_south,sound_south))//;
									+dt*MassFluxXSource(DenOld,FlxXOld,FlxYOld,0.0f,0.0f);
				FlxY[kp] = FlxYOld-(dt/dx)*(
									MassFluxYFluxX(den_east, px_east, py_east,m_east,sound_east)-
									MassFluxYFluxX(den_west, px_west, py_west,m_west,sound_west))
									-(dt/dy)*(
									MassFluxYFluxY(den_north, px_north, py_north,m_north,sound_north)-
									MassFluxYFluxY(den_south, px_south, py_south,m_south,sound_south))//;
									+dt*MassFluxYSource(DenOld,FlxXOld,FlxYOld,0.0f,0.0f);
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
	try {
		if(!isfinite(Den[pos_end]) || !isfinite(Den[pos_ini]) || !isfinite(FlxX[pos_end]) || !isfinite(FlxX[pos_ini])){
			throw "ERROR: numerical method failed to converge";
		}
		data_preview << t << "\t"
		<< Den[pos_end]  << "\t"
		<< FlxX[pos_end] << "\t"
		<< Den[pos_ini]  << "\t"
		<< FlxX[pos_ini] << "\n";
	}catch (const char* msg) {
		cerr << msg  <<"\nExiting"<< endl;
		this->CloseHdf5File();
		exit(EXIT_FAILURE);
	}
}

void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}




GrapheneFluid2D::GrapheneFluid2D(int size_nx, int size_ny, SetUpInput &input_parameters): Fluid2D(size_nx, size_ny, input_parameters){
	vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
	col_freq = input_parameters.CollisionFrequency ; // collision_frequency;
	cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency;
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2fwc=%.2f", vel_snd, vel_fer, kin_vis, col_freq,cyc_freq);
	file_infix = buffer;
}


void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

void GrapheneFluid2D::MassFluxToVelocity(){
	for(int c=0; c <= Nx * Ny - 1; c++){
		VelX[c]= FlxX[c] * pow(Den[c], -1.5f);
		VelY[c]= FlxY[c] * pow(Den[c], -1.5f);
		CurX[c] = VelX[c] * Den[c];
		CurY[c] = VelY[c] * Den[c];
	}
}



void GrapheneFluid2D::CflCondition(){ // Eventual redefinition
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );
	//dt = 2.4/(vel_snd*sqrt(25.0/(dx*dx)+16.0/(dy*dy)));
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;
	if(kin_vis>0.0f&&2.0f*kin_vis*dt > dx*dx*dy*dy/(dx*dx+dy*dy)){
		dt = 0.5f*0.25f*dx*dx/kin_vis;
	}
}	

float  GrapheneFluid2D::DensityFluxX(float n,float flx_x, __attribute__((unused)) float flx_y, __attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_1;
	f_1 = flx_x / sqrt(n);
	return f_1;
}
float  GrapheneFluid2D::DensityFluxY(float n,__attribute__((unused)) float flx_x, float flx_y, __attribute__((unused)) float mass,__attribute__((unused))  float s){
	float f_1;
	f_1 = flx_y / sqrt(n);
	return f_1;
}

float  GrapheneFluid2D::MassFluxXFluxX(float n,float flx_x,__attribute__((unused)) float flx_y,float mass, float s){
	float f_2;
	f_2 = flx_x * flx_x / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_2;
}
float  GrapheneFluid2D::MassFluxXFluxY(__attribute__((unused)) float n,float flx_x, float flx_y,float mass,__attribute__((unused)) float s){
	float f_2;
	f_2 = flx_x * flx_y / mass;
	return f_2;
}
float  GrapheneFluid2D::MassFluxYFluxX(__attribute__((unused)) float n,float flx_x, float flx_y,float mass,__attribute__((unused)) float s){
	float f_3;
	f_3 = flx_x * flx_y / mass;
	return f_3;
}
float  GrapheneFluid2D::MassFluxYFluxY(float n,__attribute__((unused)) float flx_x, float flx_y,float mass, float s){
	float f_3;
	f_3 = flx_y * flx_y / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * s * s * n * n;
	return f_3;
}

void GrapheneFluid2D::MagneticSourceSemiAnalytic(){
	float px_0,py_0,sqrtn_0;
	float wc=cyc_freq;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 * cos(wc * dt / sqrtn_0) - py_0 * sin(wc * dt / sqrtn_0);
			FlxY[kp]= px_0 * sin(wc * dt / sqrtn_0) + py_0 * cos(wc * dt / sqrtn_0);
		}
	}		
}

void GrapheneFluid2D::ViscosityFtcs() {
int north, south, east, west;
float mass_den_center, mass_den_north, mass_den_south, mass_den_east, mass_den_west;


//calculate laplacians
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		div_t divresult;
		divresult = div(kp, Nx);
		int j = divresult.quot;
		int i = divresult.rem;
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			north = i + (j + 1) * Nx;
			south = i + (j - 1) * Nx;
			east = i + 1 + j * Nx;
			west = i - 1 + j * Nx;
			mass_den_center = pow(Den[kp], 1.5f);
			mass_den_north = pow(Den[north], 1.5f);
			mass_den_south = pow(Den[south], 1.5f);
			mass_den_east = pow(Den[east], 1.5f);
			mass_den_west = pow(Den[west], 1.5f);
			lap_flxX[kp] =
					(-4.0f * FlxX[kp] / mass_den_center + FlxX[north] / mass_den_north + FlxX[south] / mass_den_south + FlxX[east] / mass_den_east +
					 FlxX[west] / mass_den_west) / (dx * dx);
			lap_flxY[kp] =
					(-4.0f * FlxY[kp] / mass_den_center + FlxY[north] / mass_den_north + FlxY[south] / mass_den_south + FlxY[east] / mass_den_east +
					 FlxY[west] / mass_den_west) / (dx * dx);
		}
	}
	//FTCS algorithm
	float old_px,old_py,sqrtn_0;
	//float	odd_vis=;
	//float	odd_vis=0.0f;
	for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
		if (kp % Nx != Nx - 1 && kp % Nx != 0) {
			old_px=FlxX[kp];
			old_py=FlxY[kp];
			sqrtn_0=sqrt(Den[kp]);
			//FlxX[kp] = old_px + dt * (kin_vis * lap_flxX[kp] - odd_vis*lap_flxY[kp]);
			//FlxY[kp] = old_py + dt * (kin_vis * lap_flxY[kp] + odd_vis*lap_flxX[kp]);
			//FlxX[kp] = old_px + dt * (kin_vis * lap_flxX[kp] -cyc_freq * old_py / sqrtn_0);
			//FlxY[kp] = old_py + dt * (kin_vis * lap_flxY[kp] +cyc_freq * old_px / sqrtn_0);
			FlxX[kp] = old_px + dt * (kin_vis * lap_flxX[kp] +col_freq*old_px*0.0f -cyc_freq * old_py / sqrtn_0);
			FlxY[kp] = old_py + dt * (kin_vis * lap_flxY[kp] +col_freq*old_px*0.0f +cyc_freq * old_px / sqrtn_0);
		}
	}
}



void GrapheneFluid2D::MagneticSourceFtcs(){
	float px_0,py_0,sqrtn_0;
	for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
		if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn_0=sqrt(Den[kp]);
			px_0=FlxX[kp];
			py_0=FlxY[kp];
			FlxX[kp]= px_0 -  dt * cyc_freq * py_0 / sqrtn_0;
			FlxY[kp]= py_0 +  dt * cyc_freq * px_0 / sqrtn_0;
		}
	}
}

float GrapheneFluid2D::DensitySource(float n, float flx_x, float flx_y, float mass, float s) {
	return 0.0f;
}
float GrapheneFluid2D::MassFluxXSource(float n, float flx_x, float flx_y, float mass, float s) {
	return -1.0f*col_freq*flx_x;
}
float GrapheneFluid2D::MassFluxYSource(float n, float flx_x, float flx_y, float mass, float s) {
	return -1.0f*col_freq*flx_y;
}



void Fluid2D::SaveSound() {
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocicity", hdf5_float, *DataspaceVelSnd);
	dataset_vel_snd.write(vel_snd_arr, hdf5_float);
	dataset_vel_snd.close();
}

void Fluid2D::SaveSnapShot() {
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );

	int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
	snapshot_step = points_per_period / snapshot_per_period;

	this->MassFluxToVelocity();
	string str_time = to_string(TimeStepCounter / snapshot_step);
	string name_dataset = "snapshot_" + str_time;

//	cout << TimeStepCounter <<"\t"<<name_dataset <<endl;


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
	dataset_vel_x.write(VelX, hdf5_float);
	dataset_vel_x.close();
	atr_step_vel_x.write( hdf5_int, &TimeStepCounter);
	atr_time_vel_x.write( hdf5_float , &currenttime);
	atr_step_vel_x.close();
	atr_time_vel_x.close();

	DataSet dataset_vel_y = GrpVelY->createDataSet(name_dataset, hdf5_float, *DataspaceVelY);
	Attribute atr_step_vel_y = dataset_vel_y.createAttribute("time step", hdf5_int, atr_dataspace);
	Attribute atr_time_vel_y = dataset_vel_y.createAttribute("time", hdf5_float, atr_dataspace);
	dataset_vel_y.write(VelY, hdf5_float);
	dataset_vel_y.close();
	atr_step_vel_y.write( hdf5_int, &TimeStepCounter);
	atr_time_vel_y.write( hdf5_float , &currenttime);
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

float Fluid2D::DensitySource(float n, float vel_x, float vel_y, float mass, float s) {
	return 0.0f;
}

float Fluid2D::MassFluxXSource(float n, float vel_x, float vel_y, float mass, float s) {
	return 0.0f;
}
float Fluid2D::MassFluxYSource(float n, float vel_x, float vel_y, float mass, float s) {
	return 0.0f;
}
