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


//TODO REVER A NECESSIDADES DESTAS FUNCOES:
//void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);
//constexpr float Gauss_Kernel(int position , float t); //
//constexpr float Gauss_Kernel_Derivative(int position , float t); //
void Record_Log_File(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float dy, float tmax);

//float Signal_Average(int n, float dt, const float * f);
float Integral_1_D(int n, float ds, const float * f);
float Integral_2_D(int n, int m, float dx, float dy, const float * f);

//void Extrema_Finding(float * vec_in, int n, float sound, float dt, float & sat, float  & tau, float & error, const std::string& extremafile);
float Imag_Freq(float sound, float fermi, float col_freq);                  //
float Phase_Vel(float sound, float fermi);
float Real_Freq(float sound, float fermi, float col_freq, int mode);  //

void Parameter_Initialization(int argc, char ** argv, int &data_save_mode, float &input_vel_snd, float &input_vel_fer, float &input_col_freq, float &input_kin_vis, float &input_cyc_freq);
void Parameter_Exceptions_Checking(const int &data_save_mode, const  float &input_vel_snd, const  float &input_vel_fer, const  float &input_col_freq, const  float &input_kin_vis, const  float &input_cyc_freq);
//-----------------------------------


/* Functions to implement the spatial variation of the sound velocity S(x) in 1D or S(x,y) in 2D
 * corresponding to a variation of substrat permitivitty or even the description of a multi gated system.
 * */
float Sound_Velocity_Anisotropy(int i, float dx, float s);
float Sound_Velocity_Anisotropy(int  i,float dx, int j,float dy, float s);
/*....................................................................................................................*/

/* Average moving filter for the smoothing of 1D simulation, suppressing the spurious oscillations inherent to the 2nd order solver*/
void Average_Filter(const float * vec_in, float * vec_out, int size , int width );


/* Base class from which the fluid and graphene fluid classes are derived
 * The setting of the dimensions and creation of the HDF5 file and its attributes is controlled from this class
 * */
class TethysBase {
	protected:
		int   Nx ;          // Simulation region (dataset) dimensions
		int   Ny ;
		int   RANK;         // Either 1 or 2 indicating the dimension of the problem
		float dx=1.0f;      // space discretization along x
		float dy=1.0f;      // space discretization along y
		float dt=1.0f;      // temporal discretization. will later be redifined by the CFL condition
		float lengX=1.0f;   // physical length along x. dx will be redifined as lengX/Nx
		float lengY=1.0f;   // physical length along y. dy will be redifined as lengY/Ny
		float vel_snd =50.0f;   // sound velocity parameter
		float kin_vis =0.0f;    // kinetic shear viscosity parameter
		float col_freq =0.0f;   // colision frequency parameter
		std::string file_infix; // base name for the output files
		float Tmax=10;          // total time of simulation
		bool HDF5fileCreated = false;   // flag to indicate the if the HDF5 file was created in the case the user choose the full output option
	public:
		TethysBase(int size_nx, int size_ny, int dimensions); // class constructor initializes Nx, Ny, RANK and file_infix
		~TethysBase();  // class destructor if the the flag HDF5fileCreated=TRUE it deletes the dataspaces and hdf5 files

		H5File* Hdf5File ;  // hdf5 file handler
		Group* GrpDat ;     // group for the simulated data: Attributes; Density; Velocity X; Velocity Y
		Group* GrpDen ;     // group for ALL Density snapshots
		Group* GrpVelX ;    // group for ALL Velocity X snapshots
		Group* GrpVelY ;    // group for ALL Velocity X snapshots
		DataSpace* DataspaceDen;    // dataspace for EACH Density snapshots
		DataSpace* DataspaceVelX;   // dataspace for EACH Velocity X snapshots
		DataSpace* DataspaceVelY;   // dataspace for EACH Velocity Y snapshots

		void SetTmax(float x);      // setter method for the total simulation time cf. GrapheneFluid2D::SetSimulationTime()
		float GetTmax() const;      // getter method for the total simulation time
		int SizeX() const;          // getter method for the number of simulation points along x
		int SizeY() const;          // getter method for the number of simulation points along y
		int Rank() const;           // getter method for the system dimensionality

		void SetVelSnd(float x);    // setter method for nominal S value
		void SetKinVis(float x);    // setter method for kinetic shear viscosity
		void SetColFreq(float x);   // setter method for collision frequency
		void SetDx(float x);        // setter method for spatial step x
		void SetDy(float x);        // setter method for spatial step y
		void SetDt(float x);        // setter method for temporal step
		void SetLengthX(float x);   // setter method for total length along x
		void SetLengthY(float x);   // setter method for total length along y

		float GetVelSnd() const;    // getter method for nominal S value
		float GetKinVis() const;    // getter method for kinetic shear viscosity
		float GetColFreq() const;   // getter method for collision frequency
		float GetDx() const;        // getter method for spatial discretization x
		float GetDy() const;        // getter method for spatial discretization y
		float GetDt() const;        // getter method for time discretization
		float GetLengthX() const;   // getter method for total length along x
		float GetLengthY() const;   // getter method for total length along y

		std::string GetInfix() const;   // getter method for file name infix
		void CreateHdf5File();          // creates the HDF5 files with the necessary structure
		void CloseHdf5File();           // closes the HDF5 file
		void WriteAttributes();          // saves the simulation attributes (either physical and simulation parameters)

		void BannerDisplay(); // launches the initial ASCII art banner
		void WelcomeScreen(float vel_fer) const; //launches screen with the relevant info
};
#endif

