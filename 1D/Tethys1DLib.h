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
		int   RANK=1;
		std::string file_infix ;
	
	public:
		TETHYSBase(int sizeN); // acho que pelo menos para já nao vai precisar de construtor ou entao ponho o banner mesmo no constrturos 
		~TETHYSBase();	

		H5File* hdf5file ; // se tirar o namespace nao esquecer usar o H5::
		Group* grp_dat ;
		Group* grp_den ;
		Group* grp_vel ;
		Group* grp_cur ;
		DataSpace* dataspace_den;
		DataSpace* dataspace_vel;
		DataSpace* dataspace_cur;
		
		std::string GetInfix();
		void virtual SetFileName();
		void CreateHDF5File();
		void BannerDisplay(void);
		void WellcomeScreen(float vel_snd, float vel_fer,float col_freq,float viscosity, float dt, float dx, float Tmax);
};  

class Fluid1D : public TETHYSBase{
	protected:	
		float dx=1.0;		
		float dt=1.0;		
		float Tmax=10;
		const float leng=1.0;
		int Nx;
		float vel_snd =50.0;
		float kin_vis =0.0;
		float * vel_snd_arr;	
		float * den_mid ;
		float * vel_mid ;
		float * grad_vel_mid;
		std::ofstream data_preview;
		
	public :
		using TETHYSBase::TETHYSBase;	 				
		float * den ;
		float * vel ;
		float * grad_vel;
		float * cur ;
		float * den_cor;
		float * vel_cor ;
		float * cur_cor ;
		
		Fluid1D(int sizeN);
		~Fluid1D();

		void SetFileName() override;
		void Smooth(int width);
		void SetVelSnd(float x);
		void SetKinVis(float x);
		float GetVelSnd();
		float GetKinVis();
		float GetDx();
		float GetDt();
		float GetTmax();
		void SetTmax(float x);
		void SetSimulationTime();
		int SizeX();
		void InitialCondRand();
		void Richtmyer();
		void SetSound();		
		virtual void CFLCondition();
		virtual float DensityFlux(float n,float v,float S);
		virtual float VelocityFlux(float n,float v,float dv,float S);
		virtual float DensitySource(float n,float v,float S);
		virtual float VelocitySource(float n,float v,float S);
		void CreateFluidFile();
		void WriteFluidFile(float t) ;
};


class GrapheneFluid1D : public Fluid1D{
	protected : 
		float vel_fer =10.0;							
		float col_freq =0.0; 		
	 	
	public : 
		using Fluid1D::Fluid1D;
	
		void SetFileName() override;
		void CFLCondition() override;
	    void BoundaryCond(int type);		
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		float DensityFlux(float n,float v,float S) override;
		float VelocityFlux(float n,float v,float dv,float S) override;
		float DensitySource(float n,float v,float S) override; 
		float VelocitySource(float n,float v,float S) override;		
		void WriteAtributes();
};



class ElectroAnalysis {
	private:
			std::ofstream data_electro;	
	public:
		void CreateElectroFile(GrapheneFluid1D& graphene);
		void WriteElectroFile(float t,GrapheneFluid1D& graphene);				
		float NetCharge(GrapheneFluid1D& graphene);
		float OhmPower(GrapheneFluid1D& graphene);
		float AverageCurrent(GrapheneFluid1D& graphene);
		float ElectricDipole(GrapheneFluid1D& graphene);
		float ElectricDipoleVariation(GrapheneFluid1D& graphene);
};

#endif

