#ifndef TESTLIB_H
#define TESTLIB_H


#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>


#include <H5Cpp.h>
using namespace std;
using namespace H5;


//REVER A NECESSIDADES DESTAS FUNCOES 
void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);
constexpr float Gauss_Kernel(int position , float t); //
constexpr float Gauss_Kernel_Derivative(int position , float t); //
void Record_Log_File(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float dy, float tmax);

float Signal_Average(int n, float dt, float * f);
float Integral_1_D(int n, float ds, float * f);

void Extrema_Finding(float * vec_in, int n, float sound, float dt, float & sat, float  & tau, float & error, std::string extremafile);
float Imag_Freq(float sound, float fermi, float col_freq);                  //
float Phase_Vel(float sound, float fermi);
float Real_Freq(float sound, float fermi, float col_freq, int mode);  //

void Parameter_Initalization(int argc, char ** argv, int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis, float &input_cyc_freq);
void Parameter_Exeptions_Checking(int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis, float &input_cyc_freq);
//-----------------------------------

float Sound_Velocity_Anisotropy(float i, float dx, float s);
float Sound_Velocity_Anisotropy(float i,float dx, float j,float dy, float s);
void Average_Filter(float * vec_in, float * vec_out, int size , int width );

class TETHYSBase {
	protected:
		int   Nx ;                    // dataset dimensions
		int   Ny ;
		int   RANK;
		std::string file_infix = "BaseFluid1D" ;
		float Tmax=10;
		bool HDF5fileCreated = false;
	public:
		TETHYSBase(int size_nx, int size_ny, int dimensions); // acho que pelo menos para já nao vai precisar de construtor ou entao ponho o banner mesmo no constrturos
		~TETHYSBase();


		H5File* hdf5file ; // se tirar o namespace nao esquecer usar o H5::
		Group* GrpDat ;
		Group* GrpDen ;
		Group* GrpVelX ;
		Group* GrpVelY ;
		DataSpace* DataspaceDen;
		DataSpace* DataspaceVelX;
		DataSpace* DataspaceVelY;

		float GetTmax();
		void SetTmax(float x);
		int SizeX();
		int SizeY();
		int Rank();
		
		std::string GetInfix();
		void CreateHDF5File();
		void CloseHDF5File();
		void BannerDisplay();
		void WellcomeScreen(float vel_snd, float vel_fer,float col_freq,float viscosity, float dt, float dx, float tmax);
		
		
};  



#endif

