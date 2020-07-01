#ifndef TETHYSLIB_1D_H
#define TETHYSLIB_1D_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;

 
class Fluid1D : public TETHYSBase{
	protected:	
		float dx=1.0;		
		float dt=1.0;		
		const float leng=1.0;
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
		
		Fluid1D(int sizeN,float VELSND, float VISCO);
		~Fluid1D();

		void Smooth(int width);
		void SetVelSnd(float x);
		void SetKinVis(float x);
		float GetVelSnd();
		float GetKinVis();
		float GetDx();
		float GetDt();
		void SetSimulationTime();
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
	
		GrapheneFluid1D(int sizeN,float VELSND, float FERMI,float VISCO,float COL);
	
		void CFLCondition() override;
//	    void BoundaryCond(int type);		
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

