#ifndef TETHYS1DLIB_H
#define TETHYS1DLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;
using namespace std;

class Fluid1D : public TETHYSBase{
	protected:
		float * vel_snd_arr;        // array for saving the (potentially varying) S(x) function
		float * den_mid ;           // number density at midpoint
		float * vel_mid ;           // velocity at midpoint
		float * grad_vel_mid;       // velocity gradient at mid point for the viscous case
		std::ofstream data_preview; // file stream for simplified .dat file output
	public :
		float * Den ;       // number density
		float * Vel ;       // fluid velocity
		float * GradVel;    // fluid velocity gradient
		float * Cur ;       // current density (density times velocity)
		float * DenCor;     // corrected i.e. smoothed quantities
		float * VelCor ;
		float * CurCor ;
		Fluid1D(int size_nx, float sound_velocity, float shear_viscosity);
		~Fluid1D();
		void Smooth(int width);     // smoothing moving average filter to obtain the "Cor" version of the quantities
		//TODO move the setters for the tethys base
		void SetVelSnd(float x);    // setter method for nominal S value
		void SetKinVis(float x);    // setter method for kinetic shear viscosity
		float GetVelSnd() const;    // getter method for nominal S value
		float GetKinVis() const;    // getter method for kinetic shear viscosity
		float GetDx() const;        // getter method for spatial discretization x
		float GetDt() const;        // getter method for time discretization

		void SetSimulationTime();   // Finds and set the appropriate simulation time that is 1) Longer than the saturation time 2) Contains enough oscillation periods in the saturated region
		void InitialCondRand();     // Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();     // Initial condition for testing and debugging
		void Richtmyer();           // Central Algorithm for solving the hyperbolic conservation law
		void SetSound();            // Applies the anisotropy to the sound velocity array
		virtual void CFLCondition();    // Calculates dx and imposes Courant–Friedrichs–Lewy condition to dt
		virtual float DensityFlux(float n,float v, __attribute__((unused)) float s);    // density equation (continuity equation) conserved flux
		virtual float VelocityFlux(float n,float v,float dv, __attribute__((unused)) float s); // velocity equation (momentum equation) conserved flux
		virtual float DensitySource( __attribute__((unused)) float n,  __attribute__((unused)) float v, __attribute__((unused)) float s); // density equation (continuity equation) source term
		virtual float VelocitySource( __attribute__((unused)) float n, __attribute__((unused)) float v, __attribute__((unused)) float s); // velocity equation (momentum equation) source term
		void CreateFluidFile();     // create and open the simplified .dat file output
		void SaveSnapShot(int time_step,int snapshot_step); // saves the all the simulated quantities on the appropriate dataspace of the HDF5 file
		void WriteFluidFile(float t) ; // writes the line of time t on the simplified .dat file output
};
class GrapheneFluid1D : public Fluid1D{
	protected : 
		float vel_fer =10.0;    // Fermi velocity of the electron fluid
	public :
		GrapheneFluid1D(int size_n, float sound_velocity, float fermi_velocity, float shear_viscosity, float collision_frequency);
		/*Override CFL condition to the case of graphene equations */
		void CFLCondition() override;
		void SetVelFer(float x);    // setter method for Fermi velocity
		float GetVelFer() const;    // getter method for Fermi velocity
		//TODO move colisio frqeucni setter and getter to base
		void SetColFreq(float x);   // setter method for colision frequency
		float GetColFreq() const;   // getter method for colision frequency
		/*Override fluxes and sources to specifics of graphene physics*/
		float DensityFlux(float n,float v,__attribute__((unused)) float s) override;
		float VelocityFlux(float n,float v,float dv,float s) override;
		float DensitySource(float n,float v,float s) override;
		float VelocitySource(float n,float v,float s) override;
};
#endif

