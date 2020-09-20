#ifndef TETHYS2DLIB_H
#define TETHYS2DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"

using namespace H5;


class Fluid2D : public TETHYSBase
{
	protected:
		float dx=1.0;
		float dy=1.0;
		float dt=1.0;
		float lengX=1.0;
		float lengY=1.0;
		float vel_snd =50.0;
		float kin_vis =0.0;
		float * vel_snd_arr;
		float * den_mid ; // 1st Aux. Grid (Nx-1)*(Ny-1)
		float * flxX_mid ;
		float * flxY_mid ;

		float * lap_flxX ; //new grids for the laplacians
		float * lap_flxY ;

		std::ofstream data_preview;
		virtual void SetFileName();

	public :
		float * Den ;
		float * VelX ;
		float * VelY ;
		float * FlxX ;
		float * FlxY ;
		float * CurX ;
		float * CurY ;
		Fluid2D(int size_nx, int size_ny, float sound_velocity, float shear_viscosity);
		~Fluid2D();
		void SetVelSnd(float x);
		void SetSound();
		float GetVelSnd() const;
		void SetKinVis(float x);
		float GetKinVis() const;
		void SetLengthX(float x);
		void SetLengthY(float x);
		float GetLengthX() const;
		float GetLengthY() const;
		float GetDx() const;
		void SetDx(float x);
		float GetDy() const;
		void SetDy(float x);
		float GetDt() const;
		void SetDt(float x);

		virtual void SetSimulationTime();

		void InitialCondRand();
		void InitialCondTest();
		void Richtmyer();	
		//void DimensionalSplittingMethod();
		virtual void CFLCondition();

		virtual float DensityFluxX(__attribute__((unused)) float n, float flx_x, __attribute__((unused)) float vel_y, __attribute__((unused)) float mass, __attribute__((unused)) float s);
		virtual float DensityFluxY(__attribute__((unused)) float n, __attribute__((unused)) float vel_x, float vel_y, __attribute__((unused)) float mass, __attribute__((unused)) float s);
		virtual float MassFluxXFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxXFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxYFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxYFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		//virtual float DensitySource(float n, float vel_x, float vel_y, float s);
		//virtual float MassFluxXSource(float n, float flx_x, float flx_y, float s);
		//virtual float MassFluxYSource(float n, float flx_x, float flx_y, float s);
		virtual void MassFluxToVelocity();
		
		void CreateFluidFile();
		void WriteFluidFile(float t) ;
};

class GrapheneFluid2D : public Fluid2D{
	protected : 
		float vel_fer =10.0f;
		float col_freq =0.0f;
		float cyc_freq =0.0f;
	public : 
		//using Fluid2D::Fluid2D;
		
		GrapheneFluid2D(int size_nx, int size_ny, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency, float cyclotron_frequency);

		
		void SetVelFer(float x);
		float GetVelFer() const;
		void SetColFreq(float x);
		float GetColFreq() const;
		float GetCycFreq() const;
		void CFLCondition() override;
		void SetSimulationTime() override;
		void MassFluxToVelocity() override;
		float DensityFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensityFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		//float DensitySource(float n, float flx_x, float flx_y, float s) override;
		//float MassFluxXSource(float n, float flx_x, float flx_y, float s)  override;
		//float MassFluxYSource(float n, float flx_x, float flx_y, float s) override;



		void MagneticSourceSemiAnalytic();
		void MagneticSourceFTCS();
		void ViscosityFTCS();
		void WriteAtributes();
		void SaveSnapShot(int time_step,int snapshot_step);
};


#endif
