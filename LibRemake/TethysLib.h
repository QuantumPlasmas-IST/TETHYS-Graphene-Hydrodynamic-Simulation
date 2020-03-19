
#ifndef TESTLIB_H
#define TESTLIB_H

//REVER A NECESSIDADES DESTAS FUNCOES 
void Autocorrelation(float * out_gamma ,float * in , int crop, int size);
void BannerDisplay(void);
void ConvolveGauss(int type, float M, float t, float * in, float * out, int size);
float GaussKernel(int position , float t); //
float GaussKernelDerivative(int position , float t); //
void RecordLogFile(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float Tmax);
float RootMeanSquare(int N, float dt, float * f);
float SignalAverage(int N, float dt, float * f);
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void WellcomeScreen(float vel_snd, float vel_fer,float col_freq, float dt, float dx, float Tmax);
void ExtremaFinding(float * vec_in, int N, float sound, float dt,float & sat, float  & tau, float & error, std::string extremafile);
float ImagFreq(float sound, float fermi, float col_freq);                  //
float PhaseVel(float sound, float fermi);
float RealFreq(float sound, float fermi, float col_freq, int mode);  //
void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile);
//-----------------------------------



float SoundVelocityAnisotropy(float i, float dx,float S);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );

class Fluid1D 
{
	protected:	
		float dx=1.0;		
		float dt=1.0;		
		float Tmax=10;
		const float leng=1.0;
		int Nx;
		float vel_snd =50.0;
		float * vel_snd_arr;	
		float * den_mid ;
		float * vel_mid ;

	public :
		float * den ;
		float * vel ;
		float * cur ;
		float * den_cor;
		float * vel_cor ;
		float * cur_cor ;
		
		explicit Fluid1D(int sizeN);
		~Fluid1D();

		void Smooth(int width);
		void SetVelSnd(float x);
		float GetVelSnd();
		float GetDx();
		void SetDx(float x);
		float GetDt();
		void SetDt(float x);
		float GetTmax();
		void SetTmax(float x);

		int SizeX();
		void InitialCondRand();
		void Richtmyer();
		void SetSound();		
		virtual void CFLCondition();
		virtual float DensityFlux(float n,float v,float S);
		virtual float VelocityFlux(float n,float v,float S);
		virtual float DensitySource(float n,float v,float S);
		virtual float VelocitySource(float n,float v,float S);
};


class GrapheneFluid1D : public Fluid1D{
	
	protected : 
		float vel_fer =10.0;							
		float col_freq =0.0; 						
	 	
	public : 
		using Fluid1D::Fluid1D;
	
		void CFLCondition();
	    void BoundaryCond(int type);		
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		float DensityFlux(float n,float v,float S);
		float VelocityFlux(float n,float v,float S);
		float DensitySource(float n,float v,float S);
		float VelocitySource(float n,float v,float S);		
};


//void HDF5SetUp(H5::H5File hdf5file, GrapheneFluid1D graph_obj )
//void HDF5SetUp( Fluid1D fluid_obj, string file_name );
//void HDF5SetUp( GrapheneFluid1D );  polimorfismo
//void HDF5SaveSnapshot( Fluid1D fluid_obj , string snap_name);

#endif
