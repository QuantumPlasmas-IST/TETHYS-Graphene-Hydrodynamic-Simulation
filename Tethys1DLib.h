#ifndef TETHYS1DLIB_H
#define TETHYS1DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;
using namespace std;

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
		float * Den ;
		float * Vel ;
		float * GradVel;
		float * Cur ;
		float * DenCor;
		float * VelCor ;
		float * CurCor ;
		Fluid1D(int size_nx, float sound_velocity, float shear_viscosity);
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
		void InitialCondTest();
		void Richtmyer();
		void SetSound();
		virtual void CFLCondition();
		virtual float DensityFlux(float n,float v, __attribute__((unused)) float s);
		virtual float VelocityFlux(float n,float v,float dv, __attribute__((unused)) float s);
		virtual float DensitySource( __attribute__((unused)) float n,  __attribute__((unused)) float v, __attribute__((unused)) float s);
		virtual float VelocitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s);
		void CreateFluidFile();
		void SaveSnapShot(int time_step,int snapshot_step);
		void WriteFluidFile(float t) ;
};
class GrapheneFluid1D : public Fluid1D{
	protected : 
		float vel_fer =10.0;
		float col_freq =0.0;
	public : 
		using Fluid1D::Fluid1D;
		GrapheneFluid1D(int size_n, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency);
		void CFLCondition() override;
		//void BoundaryCond(int type);
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		float DensityFlux(float n,float v,__attribute__((unused)) float s) override;
		float VelocityFlux(float n,float v,float dv,float s) override;
		float DensitySource(float n,float v,float s) override;
		float VelocitySource(float n,float v,float s) override;
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

