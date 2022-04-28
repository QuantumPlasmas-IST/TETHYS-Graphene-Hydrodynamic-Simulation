/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 1D general fluid
 */

#ifndef FLUID1DLIB_H
#define FLUID1DLIB_H


#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/GridLib.h"
#include "includes/StateVecLib.h"

using namespace H5;
using namespace std;

class CellHandler1D; //forward declaration



/*!
 * @brief Generalistic fluid class in one dimension, mainly for testing purposes.
 *
 * The Fluid1D class describes a regular newtonian compressible fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is constant.
 * It is maintained mainly for testing, the GrapheneFluid1D class is derived from it and overrides the necessary methods in order to describe the semi-classical electronic fluid.
 * */
class Fluid1D : public TethysBase{
	protected:

		std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 10;
		int snapshot_step = 1;


	float * vel_snd_arr;        // array for saving the (potentially varying) S(x) function



	void RichtmyerStep1();
	void RichtmyerStep2();


	void VelocityToCurrent();

	virtual float JacobianSpectralRadius( StateVec U);
	virtual float JacobianSignum(StateVec U,std::string key);


	friend class NumericalFlux;

public :

		StateVec * Umain;
		StateVec * Uaux;
		StateVec * Umid;

		float * Den ;       // number density
		float * Vel ;       // fluid velocity

		//TODO implement velocity grandient at StateVec
		float * GradVel;    // fluid velocity gradient

		float * Cur ;       // current density (density times velocity)

		explicit Fluid1D(const SetUpParameters &input_parameters);
		~Fluid1D();
		bool Snapshot() const;
//		void Smooth(int width);     ///< smoothing moving average filter to obtain the "Cor" version of the quantities
		void SetSimulationTime();   ///< Finds and set the appropriate simulation time that is 1) Longer than the saturation time 2) Contains enough oscillation periods in the saturated region
		void InitialCondRand();     ///< Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();     ///< Initial condition for testing and debugging
		void InitialCondGeneral(function<float(float)> fden, function<float(float)> fvx);

		void Richtmyer();
		void McCormack();
		void Upwind();
		void RungeKuttaTVD();
		void LaxFriedrichs();

		void CopyFields();





		void SetSound();            ///< Applies the anisotropy to the sound velocity array
		void SaveSound();
	    void SetSound(const std::function<float(float)>& func);            ///< Applies the anisotropy to the sound velocity array
		virtual void CflCondition();    ///< Calculates dx and imposes Courant–Friedrichs–Lewy condition to dt

		virtual float DensityFlux(StateVec U); ///< density equation (continuity equation) conserved flux
		virtual float VelocityFlux(StateVec U); ///< velocity equation (momentum equation) conserved flux

		virtual StateVec ConservedFlux(StateVec U);

		virtual float DensitySource(StateVec U); ///< density equation (continuity equation) source term
		virtual float VelocitySource(StateVec U); ///< velocity equation (momentum equation) source term
		virtual float DensitySource(__attribute__((unused)) float n,  __attribute__((unused)) float v, __attribute__((unused)) float s); ///< density equation (continuity equation) source term
		virtual float VelocitySource(float n, float v, float s, float d3den); ///< velocity equation (momentum equation) source term
		void CreateFluidFile();     ///< create and open the simplified .dat file output
		void SaveSnapShot(); ///< saves the all the simulated quantities on the appropriate dataspace of the HDF5 file
		void WriteFluidFile(float t) ; ///< writes the line of time t on the simplified .dat file output
		int GetSnapshotStep() const;
		int GetSnapshotFreq() const;

		void CalcDensityLaplacian(StateVec* u_vec, int size_x); ///< Calculates the laplacian of the density using finite differences method

		//float SideAverage(const float *input_array, GridPoint1D p, char side);
};




#endif

