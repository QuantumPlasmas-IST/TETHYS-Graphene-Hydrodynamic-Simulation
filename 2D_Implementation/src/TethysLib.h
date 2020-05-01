// 2D version

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



float SoundVelocityAnisotropy(float i, float dx,float j, float dy,float S);
void AverageFilter(float ** vec_in, float ** vec_out, int NpointsX , int NpointsY, int width ); // So far only working for dx = dy 

class Fluid2D 
{
	protected:	
		float dx=1.0;
		float dy=1.0;		
		float dt=1.0;		
		const float lengX=1.0;
		const float lengY;
		int Nx;
		int Ny;
		float vel_snd =50.0;
		float ** vel_snd_arr;	
		float ** den_mid ; // 1st Aux. Grid (Nx-1)*(Ny-1)
		float ** velX_mid ;
		float ** velY_mid ;
		float ** den_mid_x ; // 2nd Aux. Grid X (Nx-1)*(Ny)
		float ** velX_mid_x ;
		float ** velY_mid_x ;
		float ** den_mid_y ; // 2nd Aux. Grid Y (Nx)*(Ny-1)
		float ** velX_mid_y ;
		float ** velY_mid_y ;

	public :
		float ** den ;
		float ** velX ;
		float ** velY ;
		float ** curX ;
		float ** curY ;
		float ** den_cor ;
		float ** velX_cor ;
		float ** velY_cor ;
		float ** curX_cor ;
		float ** curY_cor ;
		
		explicit Fluid2D(int sizeNx, int sizeNy, float width);
		~Fluid2D();

		void SetVelSnd(float x);
		void SetSound();	
		float GetVelSnd();
		float GetDx();
		void SetDx(float x);
		float GetDy();
		void SetDy(float x);
		float GetDt();
		void SetDt(float x);
		int SizeX();
		int SizeY();
		void InitialCondRand();
		void Smooth(int width);
		void Richtmyer();	
		virtual void CFLCondition();
		virtual void BoundaryCond(int type);
		virtual float DensityFluxX(float n, float velX, float velY, float S);
		virtual float DensityFluxY(float n, float velX, float velY, float S);
		virtual float VelocityXFluxX(float n, float velX, float velY, float S);
		virtual float VelocityXFluxY(float n, float velX, float velY, float S);
		virtual float VelocityYFluxX(float n, float velX, float velY, float S);
		virtual float VelocityYFluxY(float n, float velX, float velY, float S);
		virtual float DensitySource(float n, float velX, float velY, float S);
		virtual float VelocityXSource(float n, float velX, float velY, float S);
		virtual float VelocityYSource(float n, float velX, float velY, float S);
};


class GrapheneFluid2D : public Fluid2D{
	
	protected : 
		float vel_fer =10.0;							
		float col_freq =0.0; 						
	 	
	public : 
		using Fluid2D::Fluid2D;
	
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		void CFLCondition();
	    void BoundaryCond(int type);		
		float DensityFluxX(float n, float velX, float velY, float S);
		float DensityFluxY(float n, float velX, float velY, float S);
		float VelocityXFluxX(float n, float velX, float velY, float S);
		float VelocityXFluxY(float n, float velX, float velY, float S);
		float VelocityYFluxX(float n, float velX, float velY, float S);
		float VelocityYFluxY(float n, float velX, float velY, float S);
		float DensitySource(float n, float velX, float velY, float S);
		float VelocityXSource(float n, float velX, float velY, float S);
		float VelocityYSource(float n, float velX, float velY, float S);
};

#endif
