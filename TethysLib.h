#ifndef TESTLIB_H
#define TESTLIB_H

#include <H5Cpp.h>

using namespace H5;

//REVER A NECESSIDADES DESTAS FUNCOES 
void ConvolveGauss(int type, float M, float t, float * in, float * out, int size);
float GaussKernel(int position , float t); //
float GaussKernelDerivative(int position , float t); //
void RecordLogFile(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float Tmax);

float SignalAverage(int N, float dt, float * f);
float Integral1D(int N, float ds, float * f);

void ExtremaFinding(float * vec_in, int N, float sound, float dt,float & sat, float  & tau, float & error, std::string extremafile);
float ImagFreq(float sound, float fermi, float col_freq);                  //
float PhaseVel(float sound, float fermi);
float RealFreq(float sound, float fermi, float col_freq, int mode);  //

//-----------------------------------

float SoundVelocityAnisotropy(float i, float dx,float S);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );

class TETHYSBase {
	protected:
		int   Nx ;                    // dataset dimensions
		int   Ny ;           
		int   RANK;
		std::string file_infix = "BaseFluid1D" ;
		float Tmax=10;
	public:
		TETHYSBase(int sizeNX,int sizeNY,int dimensions); // acho que pelo menos para já nao vai precisar de construtor ou entao ponho o banner mesmo no constrturos 
//		~TETHYSBase();	

		H5File* hdf5file ; // se tirar o namespace nao esquecer usar o H5::
		Group* grp_dat ;
		Group* grp_den ;
		Group* grp_velX ;
		Group* grp_velY ;
		DataSpace* dataspace_den;
		DataSpace* dataspace_velX;
		DataSpace* dataspace_velY;
		
		float GetTmax();
		void SetTmax(float x);
		int SizeX();
		int SizeY();
		int Rank();
		
		std::string GetInfix();
		void CreateHDF5File();
		void BannerDisplay(void);
		void WellcomeScreen(float vel_snd, float vel_fer,float col_freq,float viscosity, float dt, float dx, float Tmax);
};  



#endif

