#ifndef TETHYSLIB_2D_H
#define TETHYSLIB_2D_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;

class Fluid2D : public TETHYSBase
{
	protected:
		float dx=1.0;
		float dy=1.0;
		float dt=1.0;
		const float lengX=1.0;
		const float lengY=1.0;
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
		explicit Fluid2D(int size_nx, int size_ny, float sound_velocity, float shear_viscosity);
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

		virtual void SetSimulationTime();

		void InitialCondRand();
		void InitialCondTest();
		void Richtmyer();	
		//void DimensionalSplittingMethod();
		virtual void CFLCondition();

		virtual float DensityFluxX(float n, float flx_x, float vel_y, float mass, float s);
		virtual float DensityFluxY(float n, float vel_x, float vel_y,float mass, float s);
		virtual float MassFluxXFluxX(float n, float flx_x, float flx_y,float mass, float s);
		virtual float MassFluxXFluxY(float n, float flx_x, float flx_y,float mass, float s);
		virtual float MassFluxYFluxX(float n, float flx_x, float flx_y,float mass, float s);
		virtual float MassFluxYFluxY(float n, float flx_x, float flx_y,float mass, float s);
		virtual float DensitySource(float n, float vel_x, float vel_y, float s);
		virtual float MassFluxXSource(float n, float flx_x, float flx_y, float s);
		virtual float MassFluxYSource(float n, float flx_x, float flx_y, float s);
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
		
		GrapheneFluid2D(int size_nx, int size_ny, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency);

		
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		void CFLCondition() override;
		void SetSimulationTime() override;
		void MassFluxToVelocity() override;
		float DensityFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensityFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensitySource(float n, float flx_x, float flx_y, float s) override;
		float MassFluxXSource(float n, float flx_x, float flx_y, float s)  override;
		float MassFluxYSource(float n, float flx_x, float flx_y, float s) override;
		
		void MagneticSource();
		void SourceFTCS();
		void ViscosityFTCS();
		void WriteAtributes();
};
#endif
