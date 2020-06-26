// 2D version

#ifndef TETHYS2DLIB_H
#define TETHYS2DLIB_H

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
void AverageFilter(float * vec_in, float * vec_out, int NpointsX , int NpointsY, int width ); // So far only working for dx = dy 

class Fluid2D 
{
	protected:	
		float dx=1.0;
		float dy=1.0;		
		float dt=1.0;		
		const float lengX=1.0;
		const float lengY=1.0;
		int Nx;
		int Ny;
		float vel_snd =50.0;
		float kin_vis =0.0;
		float * vel_snd_arr;	
		float Tmax=10;		
		float * den_mid ; // 1st Aux. Grid (Nx-1)*(Ny-1)
		float * flxX_mid ;
		float * flxY_mid ;		
		float * den_mid_x ; // 2nd Aux. Grid X (Nx-1)*(Ny)
		float * flxX_mid_x ;
		float * flxY_mid_x ;
		float * den_mid_y ; // 2nd Aux. Grid Y (Nx)*(Ny-1)
		float * flxX_mid_y ;
		float * flxY_mid_y ;

		void RunBorderY(int j);
		void RunBorderX(int i);

		std::string file_infix ;
		std::ofstream data_preview;
		virtual void SetFileName();

	public :	
		float * den ;
		float * velX ;
		float * velY ;
		float * flxX ;
		float * flxY ;
		float * curX ;
		float * curY ;
		
		explicit Fluid2D(int sizeNx, int sizeNy, float width);
		~Fluid2D();
		void SetVelSnd(float x);
		void SetSound();	
		float GetVelSnd();
		void SetKinVis(float x);
		float GetKinVis();
		float GetDx();
		void SetDx(float x);
		float GetDy();
		void SetDy(float x);
		float GetDt();
		void SetDt(float x);
		float GetTmax();
		void SetTmax(float x);

		void SetSimulationTime();

		int SizeX();
		int SizeY();
		void InitialCondRand();
		void InitialCondTEST();
//		void Smooth(int width); NOT YET IMPLEMENTED
		void Richtmyer();	
		virtual void CFLCondition();
		virtual void BoundaryCond(int type);
		
		//virtual void BCOpenAllX(int i);
		//virtual void BCOpenAllY(int j);
		//void BCOpenDensityX(int i);
		//void BCOpenDensityY(int j); 
		//virtual void BCOpenVelocityX(int i);
		//virtual void BCOpenVelocityY(int j); 
		//void BCPeriodicX(int i);
		//void BCPeriodicY(int j);
		//void BCDirichletDensityX(int i,float density);
		//void BCDirichletDensityY(int j,float density);
		//virtual void BCNormalVelocityX(int i,float vel_nor);
		//virtual void BCNormalVelocityY(int j,float vel_nor);
		//virtual void BCTangentVelocityX(int i,float vel_nor);
		//virtual void BCTangentVelocityY(int j,float vel_nor);
	
		virtual float DensityFluxX(float n, float velX, float velY, float S);
		virtual float DensityFluxY(float n, float velX, float velY, float S);
		virtual float MassFluxXFluxX(float n, float flxX, float flxY, float S);
		virtual float MassFluxXFluxY(float n, float flxX, float flxY, float S);
		virtual float MassFluxYFluxX(float n, float flxX, float flxY, float S);
		virtual float MassFluxYFluxY(float n, float flxX, float flxY, float S);
		virtual float DensitySource(float n, float flxX, float flxY, float S);
		virtual float MassFluxXSource(float n, float flxX, float flxY, float S);
		virtual float MassFluxYSource(float n, float flxX, float flxY, float S);
		virtual void MassFluxToVelocity();
		
		void CreateFluidFile();
		void WriteFluidFile(float t) ;
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
		void CFLCondition() override;
	    void BoundaryCond(int type) override ;		
	    
	    void SetSimulationTime();
	    
		//void BCOpenAllX(int i) override;
		//void BCOpenAllY(int j) override; 
		//void BCOpenVelocityX(int i) override;
		//void BCOpenVelocityY(int j) override; 
		//void BCConstantCurrentX(int i,float current);
		//void BCConstantCurrentY(int j,float current);
	    
	    void MassFluxToVelocity() override;	    
		float DensityFluxX(float n, float flxX, float flxY, float S) override;
		float DensityFluxY(float n, float flxX, float flxY, float S) override;
		float MassFluxXFluxX(float n, float flxX, float flxY, float S) override;
		float MassFluxXFluxY(float n, float flxX, float flxY, float S) override;
		float MassFluxYFluxX(float n, float flxX, float flxY, float S) override;
		float MassFluxYFluxY(float n, float flxX, float flxY, float S) override;
		float DensitySource(float n, float flxX, float flxY, float S) override;
		float MassFluxXSource(float n, float flxX, float flxY, float S)  override;
		float MassFluxYSource(float n, float flxX, float flxY, float S) override;
};
#endif
