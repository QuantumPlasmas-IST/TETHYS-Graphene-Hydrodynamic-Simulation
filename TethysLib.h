#ifndef TETHYSLIB_H
#define TETHYSLIB_H


#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <random>
#include <exception>

#include <H5Cpp.h>
using namespace std;
using namespace H5;


//REVER A NECESSIDADES DESTAS FUNCOES 
void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);
constexpr float Gauss_Kernel(int position , float t); //
constexpr float Gauss_Kernel_Derivative(int position , float t); //
void Record_Log_File(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float dy, float tmax);

float Signal_Average(int n, float dt, const float * f);
float Integral_1_D(int n, float ds, const float * f);
float Integral_2_D(int n, int m, float dx, float dy, const float * f);

void Extrema_Finding(float * vec_in, int n, float sound, float dt, float & sat, float  & tau, float & error, const std::string& extremafile);
float Imag_Freq(float sound, float fermi, float col_freq);                  //
float Phase_Vel(float sound, float fermi);
float Real_Freq(float sound, float fermi, float col_freq, int mode);  //

void Parameter_Initialization(int argc, char ** argv, int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis, float &input_cyc_freq);
void Parameter_Exeptions_Checking(const int &data_save_mode,const  float &input_vel_snd,const  float &input_vel_fer,const  float &input_col_freq,const  float &input_kin_vis,const  float &input_cyc_freq);
//-----------------------------------

float Sound_Velocity_Anisotropy(int i, float dx, float s);
float Sound_Velocity_Anisotropy(int  i,float dx, int j,float dy, float s);
void Average_Filter(const float * vec_in, float * vec_out, int size , int width );

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

		float GetTmax() const;
		void SetTmax(float x);
		int SizeX() const;
		int SizeY() const;
		int Rank() const;
		
		std::string GetInfix() const ;
		void CreateHDF5File();
		void CloseHDF5File();
		void BannerDisplay();
		void WellcomeScreen(float vel_snd, float vel_fer,float col_freq,float viscosity, float dt, float dx,float dy, float tmax);
};

class TethysException : public exception
{
public:
	virtual const char* what() const noexcept { return "ERROR: numerical method failed to converge";}
};

class PhysicalException : public exception
{
public:
	virtual const char* what() const noexcept { return "ERROR: Unphysical situation reached";}
};

class ParameterException : public exception
{
public:
	virtual const char* what() const noexcept { return "ERROR: Unphysical parameter";}
};


/*
 *
class My_Exception : public std::exception
{
public:
virtual char const * what() const noexcept { return "Something bad happend."; }
};

 try {
    ComplexOperationThatCouldFailABunchOfWays();
} catch (std::exception& e) {
    cerr << e.what() << endl;
}


*/
#endif

