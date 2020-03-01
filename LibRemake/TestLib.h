#ifndef TESTLIB_H
#define TESTLIB_H

float SoundVelocityAnisotropy(float i, float dx,float S);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );

class Fluid1D 
{
	protected:	
		float dx=1.0;		
		float dt=1.0;		
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

#endif
