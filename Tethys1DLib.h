#ifndef TETHYS1DLIB_H
#define TETHYS1DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;
using namespace std;

class Fluid1D : public TETHYSBase{
	protected:
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
		float GetVelSnd() const;
		float GetKinVis() const;
		float GetDx() const;
		float GetDt() const;
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
	public :
		GrapheneFluid1D(int size_n, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency);
		void CFLCondition() override;
		void SetVelFer(float x);
		float GetVelFer() const;
		void SetColFreq(float x);
		float GetColFreq() const;
		float DensityFlux(float n,float v,__attribute__((unused)) float s) override;
		float VelocityFlux(float n,float v,float dv,float s) override;
		float DensitySource(float n,float v,float s) override;
		float VelocitySource(float n,float v,float s) override;
};
#endif

