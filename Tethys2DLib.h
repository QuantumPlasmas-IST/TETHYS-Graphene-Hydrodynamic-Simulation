#ifndef TETHYS2DLIB_H
#define TETHYS2DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"

using namespace H5;


class Fluid2D : public TethysBase
{
	protected:
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
		void SetSound();


		virtual void SetSimulationTime();

		void InitialCondRand();
		void InitialCondTest();
		void Richtmyer();
		virtual void CflCondition();

		virtual float DensityFluxX(__attribute__((unused)) float n, float flx_x, __attribute__((unused)) float vel_y, __attribute__((unused)) float mass, __attribute__((unused)) float s);
		virtual float DensityFluxY(__attribute__((unused)) float n, __attribute__((unused)) float vel_x, float vel_y, __attribute__((unused)) float mass, __attribute__((unused)) float s);
		virtual float MassFluxXFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxXFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxYFluxX(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);
		virtual float MassFluxYFluxY(float n, float flx_x, float flx_y,__attribute__((unused)) float mass, float s);

		virtual void MassFluxToVelocity(); // Converts the mass flux density p=mnv to velocity
		
		void CreateFluidFile();
		void WriteFluidFile(float t) ;
};

class GrapheneFluid2D : public Fluid2D{
	protected : 
		float vel_fer =10.0f;
		float cyc_freq =0.0f;
	public :
		GrapheneFluid2D(int size_nx, int size_ny, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency, float cyclotron_frequency);

		void SetVelFer(float x);
		float GetVelFer() const;

		float GetCycFreq() const;
		void CflCondition() override;
		void SetSimulationTime() override;
		void MassFluxToVelocity() override; // Converts the mass density flux back to velocity, in graphene  v = p n^{-3/2}
		/*Override fluxes and sources to specifics of graphene physics*/
		float DensityFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float DensityFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxXFluxY(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxX(float n, float flx_x, float flx_y,float mass, float s) override;
		float MassFluxYFluxY(float n, float flx_x, float flx_y,float mass, float s) override;

		void MagneticSourceSemiAnalytic();
		void MagneticSourceFtcs();
		void ViscosityFtcs();
		void SaveSnapShot(int time_step,int snapshot_step);
};


#endif
