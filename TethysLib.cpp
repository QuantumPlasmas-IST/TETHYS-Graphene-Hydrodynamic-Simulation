#include "TethysLib.h"


using namespace H5;
using namespace std;

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846f
#endif

#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905f
#endif

void Parameter_Exeptions_Checking(int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis, float &input_cyc_freq){
	if(input_vel_snd<=0.0f){
		throw "Unphysical Sound Velocity";
	}
	if(input_vel_fer<=0.0f){
		throw "Unphysical Fermi Velocity";
	}
	if(input_kin_vis<0.0f){
		throw "Unphysical Shear Viscosity";
	}
	if(input_col_freq<0.0f){
		throw "Unphysical Collision Frequency";
	}
	if(input_cyc_freq<0.0f){
		throw "Unphysical Cyclotron Frequency";
	}
	if( data_save_mode != 0 && data_save_mode != 1  ) {
		throw "Unknown save mode option";
	}
}

void Parameter_Initalization(int argc, char ** argv, int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis,float &input_cyc_freq){
	if(argc==7){
		try {
			input_vel_snd = static_cast<float>(atof(argv[1]));
			input_vel_fer = static_cast<float>(atof(argv[2]));
			input_col_freq = static_cast<float>(atof(argv[3]));
			input_kin_vis = static_cast<float>(atof(argv[4]));
			input_cyc_freq = static_cast<float>(atof(argv[5]));
			data_save_mode = atoi(argv[6]);    // full data or light save option
			Parameter_Exeptions_Checking(data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis,input_cyc_freq);
		}catch (const char* msg) {
			cerr << msg << endl;
			exit(EXIT_FAILURE);
		}
	}
	else{
		try {
			cout << "Define S value: "; // throw exceptions if the velocities or frequency are negative or if S<Vf
			cin >> input_vel_snd;
			cout << "Define vF value: ";
			cin >> input_vel_fer;
			cout << "Define kinetic viscosity: ";
			cin >> input_kin_vis;
			cout << "Define collision frequency: ";
			cin >> input_col_freq;
			cout << "Define cyclotron frequency: ";
			cin >> input_cyc_freq;
			cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
			cin >> data_save_mode;
			Parameter_Exeptions_Checking(data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis, input_cyc_freq);
		}catch (const char* msg) {
			cerr << msg << endl;
			exit(EXIT_FAILURE);
		}
	}
}



float Sound_Velocity_Anisotropy(float i, float dx, float s){
	return s;
}
float Sound_Velocity_Anisotropy(float i,float dx, float j,float dy, float s){
	return s;
}


void Average_Filter(float * vec_in, float * vec_out, int size , int width ){
	for ( int i = 0; i < size; i++ ){
		if(i>=width &&i<=size-1-width){
			for(int k = i-width; k <= i+width;k++){
				vec_out[i] += vec_in[k]; 
			}
			vec_out[i] = vec_out[i]/(2.0f*(float)width+1.0f);
			}
		else{
			vec_out[i] = vec_in[i] ;
		}
	}
}





/*....................................................................*/
/*........ General Functions .........................................*/
/*....................................................................*/
void TETHYSBase::BannerDisplay(){
cout<<"\n" ;
	cout<<"╔═════════════════════════════════════════════════════════════════════════╗\n";
	cout<<"║\033[2m  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆▆▆▆▆▆  ▆▆▆▆▆▆▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▆▆ ▆▆▆▖   ▗▆▆▆ ▗▆▆▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m  █▘  ▐█▌  ▝█  ▐█▌    ▝█  █▘  ▐█▌  ▝█  ▐█▌   ▐█▌     █▌ ▐█   ▐█▌     ▝█  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌▆▆▆█        ▐█▌      ▐█▌▆▆▆▐█▌      ▐█▌     ▝██▆▆▆▆▆▖  \033[0m║\n";
	cout<<"║\033[2m      ▐█▌      ▐█▌    ▗▉      ▐█▌      ▐█▌   ▐█▌      ▐█▌    ▗       ██  \033[0m║\n";
	cout<<"║\033[2m     ▆███▆    ▆███▆▆▆██▉     ▆███▆    ▆███▆ ▆███▆    ▆███▆   ▐█▆▆▆▆▆██▘  \033[0m║\n";
	cout<<"║                                                                         ║\n";
	cout<<"║ \033[1mTwo-dimensional Emitter of THz, Hydrodynamic Simulation.  Version 2.1.0\033[0m ║\n";
	cout<<"╚═════════════════════════════════════════════════════════════════════════╝\n";
}

void TETHYSBase::WellcomeScreen(float vel_snd, float vel_fer, float col_freq,float viscosity, float dt,float dx, float tmax){
	cout << "\nFermi velocity\t\033[1mvF\t"<< vel_fer <<" v\342\202\200\033[0m\n";
	if (Phase_Vel(vel_snd, vel_fer) < vel_fer){
		cout << "Phase velocity\t\033[1mS'\t" << Phase_Vel(vel_snd, vel_fer) << " v\342\202\200\033[0m  \033[1;5;7;31m WARNING plasmon in damping region \033[0m" << endl;
	}else{
		cout << "Phase velocity\t\033[1mS'\t" << Phase_Vel(vel_snd, vel_fer) << " v\342\202\200\033[0m\n";
	}
	cout << "Viscosity \t\033[1m\316\267\t"<< viscosity <<"\033[0m\n";
	if (viscosity !=0.0){ cout << "Reynolds n. \t\033[1mRe\t"<< 1.0/viscosity <<"\033[0m\n";}
	cout << "Collision \t\033[1m\316\275\t"<< col_freq <<" v\342\202\200/L\n\033[0m\n";
	cout << "Theoretical frequency \033[1m\317\211=\317\211'+i\317\211''\033[0m\n";
	cout << "\033[1m\317\211'\t" << Real_Freq(vel_snd, vel_fer, col_freq, 1) << " v\342\202\200/L\t2\317\200/\317\211'\t" << 2.0 * MAT_PI /
	                                                                                                                         Real_Freq(vel_snd, vel_fer, col_freq, 1) << " L/v\342\202\200\033[0m\n";
	cout << "\033[1m\317\211''\t" << Imag_Freq(vel_snd, vel_fer, col_freq) << " v\342\202\200/L\t2\317\200/\317\211''\t" << 2.0 * MAT_PI /
	                                                                                                                        Imag_Freq(vel_snd, vel_fer, col_freq) << " L/v\342\202\200\033[0m\n";
	cout << "Determined maximum simulated time\t\033[1m\nT\342\202\230\342\202\220\342\202\223\t" << tmax << " L/v\342\202\200\033[0m\n";
	cout <<"Discretisation\n";
	cout <<"\033[1m\316\224t\t"<<dt<<" L/v\342\202\200\t\316\224x\t"<<dx<<" L\033[0m\n"<<endl;
}

std::string TETHYSBase::GetInfix(){return file_infix;}


float TETHYSBase::GetTmax(){return Tmax;}
void TETHYSBase::SetTmax(float x){ Tmax=x;}
int TETHYSBase::Rank(){ return RANK; }
int TETHYSBase::SizeX(){ return Nx; }
int TETHYSBase::SizeY(){ return Ny; }


void TETHYSBase::CloseHDF5File(){
	GrpDat->close();
	GrpDen->close();
	GrpVelX->close();
	if (RANK == 2) {
		GrpVelY->close();
	}
	hdf5file->close();
}


TETHYSBase::~TETHYSBase(){
	if(HDF5fileCreated) {
		delete GrpDat;
		delete GrpDen;
		delete GrpVelX;
		delete DataspaceDen;
		delete DataspaceVelX;
		if (RANK == 2) {
			delete GrpVelY;
			delete DataspaceVelY;
		}
		delete hdf5file;
	}
}









TETHYSBase::TETHYSBase(int size_nx, int size_ny, int dimension){
	Nx = size_nx;
	Ny = size_ny;
	RANK=dimension;
	const FloatType hdf5_float(PredType::NATIVE_FLOAT);
	const IntType hdf5_int(PredType::NATIVE_INT);
	char buffer [50];
	sprintf (buffer, "Fluido1D_Nx=%d", Nx);
	file_infix = buffer;
}

void TETHYSBase::CreateHDF5File(){
	std::string hdf5name;
	HDF5fileCreated=true;
	if(RANK==1){
		hdf5name = "hdf5_1D_" + this->GetInfix() + ".h5" ;
	}
	if(RANK==2){
		hdf5name = "hdf5_2D_" + this->GetInfix() + ".h5" ;
	}
	H5std_string  FILE_NAME( hdf5name );
	hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

	GrpDat = new Group(hdf5file->createGroup("/Data" ));
	GrpDen = new Group(hdf5file->createGroup("/Data/Density" ));
	GrpVelX = new Group(hdf5file->createGroup("/Data/VelocityX" ));
	if(RANK==2) {
		GrpVelY = new Group(hdf5file->createGroup("/Data/VelocityY"));
	}
	if(RANK==1){		
		hsize_t dimsf[1];// dataset dimensions
		dimsf[0] = static_cast<hsize_t>(Nx);
		DataspaceDen = new DataSpace(RANK, dimsf );
		DataspaceVelX = new DataSpace(RANK, dimsf );
	}
	if(RANK==2){
		hsize_t dimsf[2];
		dimsf[0] = static_cast<hsize_t>(Ny);
		dimsf[1] = static_cast<hsize_t>(Nx);  //troquei !
		DataspaceDen = new DataSpace(RANK, dimsf );
		DataspaceVelX = new DataSpace(RANK, dimsf );
		DataspaceVelY = new DataSpace(RANK, dimsf );
	}
}


void Record_Log_File(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float dy, float tmax){
	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	char buffer [80];
	strftime (buffer,80,"%F %H:%M:%S\n",time_info);
	logfile << "\n#Simulation @ " << buffer ;
	logfile << "#parameters:\n";
	logfile << "#vel_snd \t vel_fer \t col_freq  \t w' \t w'' \n";
	logfile << vel_snd << "\t" << vel_fer << "\t" << col_freq << "\t" << Real_Freq(vel_snd, vel_fer, col_freq, 1) << "\t" << Imag_Freq(
			vel_snd, vel_fer, col_freq) << "\n";
	logfile << "#discretisation:\n";
	logfile << "#dt\tdx\ttmax\ttime steps\tspace points\n";
	logfile << dt << "\t" << dx << "\t" << dy << "\t" << tmax << "\t" << (int) tmax / dt << "\t" << (int) 1 / dx << endl;
}

float Integral_1_D(int n, float ds, float * f){
	float itg=0.0;
	for(int j=1; j < n / 2; j++){
		itg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	itg = itg*ds/3.0f;
	return itg;	
}

float Signal_Average(int n, float dt, float * f){
	float avg=0.0;
	for(int j=1; j < n / 2; j++){
		avg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	avg = avg*dt/3.0f;
	avg = avg/(n * dt);
	return avg;
}

constexpr float Gauss_Kernel(int position , float t){
	return exp(-0.5f*position*position/t)/(sqrt(2.0f*MAT_PI*t));
}


constexpr float Gauss_Kernel_Derivative(int position , float t){
	return (-position*exp(-0.5f*position*position/t)/t)/(sqrt(2.0f*MAT_PI*t));
}


void Convolve_Gauss(int type, int m, float t, float * in, float * out, int size){
	if(type==0){
		for(int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-m; k <= m; k++){
					out[i] += in[i-k] * Gauss_Kernel(k, t);
				}
			}
		}
	}
	if(type==1){
		for(int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-m; k <= m; k++){
					out[i] += in[i-k] * Gauss_Kernel_Derivative(k, t);
				}
			}
			//out[i] = out[i] * size;
		}
	}
}

float Phase_Vel(float sound, float fermi){
	float vel_phs = sqrt(sound*sound+0.5f*fermi*fermi + 0.0625f );
	return vel_phs ;
}

float Real_Freq(float sound, float fermi, float col_freq, int mode){
	float vel_phs = Phase_Vel(sound, fermi);
	float vel_phs_sqr = vel_phs*vel_phs ;
	if (1 < vel_phs ){
		mode = 2*mode-1;
	}
	else{
		mode = 2*mode;
		}
	return fabs(vel_phs_sqr - 0.5625f ) * MAT_PI * mode / (2.0f * vel_phs );
}
float Imag_Freq(float sound, float fermi, float col_freq){
	float vel_phs = Phase_Vel(sound, fermi);
	float vel_phs_sqr = vel_phs*vel_phs ;
	return (vel_phs_sqr - 0.5625f ) * log(fabs( (vel_phs+0.75f)/(vel_phs-0.75f) )) / (2.0f * vel_phs ) - col_freq*(1.0f-0.125f/vel_phs);
}


void Extrema_Finding(float *vec_in, int n, float sound, float dt, float & sat, float & tau , float & error, std::string extremafile){
//TODO review the entire function
	ofstream data_extrema;
	data_extrema.open(extremafile);
	data_extrema << "#pos_max" << "\t" << "Max" <<"\t"<< "pos_min" <<"\t"<< "Min"<<endl;
	int w = static_cast<int>(floor(1.2f * 2.0f * MAT_PI / (Real_Freq(sound, 1.0f, 1.0f, 1) * dt)));
	int k = 0;
	int m = static_cast<int>(ceil(0.5f * dt * n * Real_Freq(sound, 1.0f, 1.0f, 1) / MAT_PI));
	float *vec_max;
	vec_max =(float*) calloc (static_cast<size_t>(m), sizeof(float));
	float *vec_pos;
	vec_pos =(float*) calloc (static_cast<size_t>(m), sizeof(float));
	for(int shift=0; shift < n - w ; shift += w){
		float maximum =  *max_element( vec_in + shift , vec_in + shift + w );
		int pos_max = static_cast<int>(max_element(vec_in + shift, vec_in + shift + w) - vec_in);
		float minimum =  *min_element( vec_in + shift , vec_in + shift + w );
		int pos_min = static_cast<int>(min_element(vec_in + shift, vec_in + shift + w) - vec_in);
		data_extrema << pos_max*dt << "\t" << maximum <<"\t"<< pos_min*dt <<"\t"<< minimum  <<endl;
		vec_max[k] = maximum;
		vec_pos[k] = pos_max*dt;
		k++;
	}
	sat =  *max_element( vec_max , vec_max + m );
	for(int i=1; i < m - 1; i++){
		if( vec_max[i] > 0.99f*sat ){
			tau  = 	vec_pos[i];
			error = 0.5f*(vec_pos[i+1]-vec_pos[i]);
			break;
		}
	}
	data_extrema.close();
}
