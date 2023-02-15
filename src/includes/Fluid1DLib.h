/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
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
#include "includes/StateVec1DLib.h"

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



		float * vel_snd_arr;        // array for saving the (potentially varying) S(x) function
	    float * vel_snd_arr_mid;        // array for saving the (potentially varying) S(x) function
		float * den_mid ;           // number density at midpoint
		float * vel_mid ;           // velocity at midpoint
		float * grad_vel_mid;       // velocity gradient at mid point for the viscous case
		std::ofstream data_preview; // file stream for simplified .dat file output
		int snapshot_per_period = 10;
		int snapshot_step = 1;

	float * lap_den_mid;
	float * lap_den;
	float * d3_den_mid;
	float * d3_den;

	//auxiliary pointers
	float *ptr_snd;
	float *ptr_den;
	float *ptr_vel;
	float *ptr_dendx;
	float *ptr_veldx;
	//float *ptr_tmp;
	float *ptr_lap_den;


//	virtual void BohmPotencial(string grid);
//	virtual void BohmSource(string grid);

//	int HopscotchFunction(const gsl_vector *x, gsl_vector *f);
//	static int gslwrapperHopscotchFunction(const gsl_vector *x, void *p, gsl_vector *f);

	void RichtmyerStep1();
	void RichtmyerStep2();
	void VelocityToCurrent();

//	gsl_matrix * BTCSmatrix ;
//	int permutation_index_s;
//	gsl_permutation * permutation_matrix ;

	virtual float JacobianSpectralRadius(StateVec1D U);
	friend class NumericalFlux;

public :

	StateVec1D * Umain;
	StateVec1D * Uaux;

		float * Den ;       // number density
		float * Vel ;       // fluid velocity
		float * GradVel;    // fluid velocity gradient
		float * Cur ;       // current density (density times velocity)
		float * DenCor;     // corrected i.e. smoothed quantities
		float * VelCor ;
		float * CurCor ;
		explicit Fluid1D(const SetUpParameters &input_parameters);
		~Fluid1D();
		bool Snapshot() const;
		void Smooth(int width);     ///< smoothing moving average filter to obtain the "Cor" version of the quantities
		void SetSimulationTime();   ///< Finds and set the appropriate simulation time that is 1) Longer than the saturation time 2) Contains enough oscillation periods in the saturated region
		void InitialCondRand();     ///< Initial condition, zero velocity and constant density with 0.5% white noise
		void InitialCondTest();     ///< Initial condition for testing and debugging
		void Richtmyer(); ///< Central Algorithm for solving the hyperbolic conservation law
	void McCormack();
	void Upwind();
	void LaxFriedrichs();



void RungeKuttaTVD();


	void BohmOperator(float bohm);
		void SetSound();            ///< Applies the anisotropy to the sound velocity array
	    void SetSound(std::function<float(float)> func);            ///< Applies the anisotropy to the sound velocity array
		virtual void CflCondition();    ///< Calculates dx and imposes Courant–Friedrichs–Lewy condition to dt
	//	virtual float DensityFlux(float n,float v, __attribute__((unused)) float s);    ///< density equation (continuity equation) conserved flux
	//	virtual float VelocityFlux(float n, float v, float dv, float s, float d2n); ///< velocity equation (momentum equation) conserved flux

	virtual float DensityFlux(GridPoint1D p, char side);    ///< density equation (continuity equation) conserved flux
	virtual float DensityFlux(StateVec1D U);
	virtual float VelocityFlux(GridPoint1D p, char side); ///< velocity equation (momentum equation) conserved flux
	virtual float VelocityFlux(StateVec1D U);
	virtual StateVec1D ConservedFlux(StateVec1D U);

//	void SetBTCSmatrix(float beta);

		virtual float DensitySource(__attribute__((unused)) float n,  __attribute__((unused)) float v, __attribute__((unused)) float s); ///< density equation (continuity equation) source term
		virtual float VelocitySource(float n, float v, float s, float d3den); ///< velocity equation (momentum equation) source term
		void CreateFluidFile();     ///< create and open the simplified .dat file output
		void SaveSnapShot(); ///< saves the all the simulated quantities on the appropriate dataspace of the HDF5 file
		void WriteFluidFile(float t) ; ///< writes the line of time t on the simplified .dat file output
		int GetSnapshotStep() const;
		int GetSnapshotFreq() const;

		void ChooseGridPointers(const string &grid);
		float SideAverage(const float *input_array, GridPoint1D p, char side);
};




#endif

