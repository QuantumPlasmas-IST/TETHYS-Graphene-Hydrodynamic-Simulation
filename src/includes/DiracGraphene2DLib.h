/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 2D Dirac-Fermi fluid
 */

#ifndef DIRACGRAPHENE2DLIBNOT_H
#define DIRACGRAPHENE2DLIBNOT_H

#include "includes/Fluid2DLib.h"



/*!
 * @brief Graphene electronic fluid class in two dimensions.
 *
 * The DiracGraphene2D class describes a  fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is *not* constant.
 * It overrides class Fluid2D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class DiracGraphene2D : public Fluid2D{
protected:
	/*!
	 * @brief Converts the number density to efective mass
	 *
	 * Since the mass of the fluid element is not a constant the in the graphene electronic fluid, one needs to perform the transformation
	   @f[ m^\star = n^{3/2} @f]
	 * */
	float DensityToMass(float density) override;

	StateVec2D *ptr_StateVecHole = nullptr;


//	void ForwardTimeOperator() override; ///< Time evolution for the FTCS method employed for the parabolic operators.
	virtual void ChooseGridPointers(const string &grid) override;

	float vel_therm = 10.0f; //new constant - pressure term
	float A = 0.1f; //new constant - source function, equilibrium relaxation
	float B = 1.0f; //new constant - source function, electron-hole creation

public :

	StateVec2D * HoleUmain;
	StateVec2D * HoleUmid;

	Group* GrpHDen ;     ///< group for ALL Hole Density snapshots
	Group* GrpHVelX ;    ///< group for ALL Hole Velocity X snapshots
	Group* GrpHVelY ;    ///< group for ALL Holes Velocity X snapshots
	DataSpace* DataspaceHDen;    ///< dataspace for EACH Hole Density snapshots
	DataSpace* DataspaceHVelX;   ///< dataspace for EACH Hole Velocity X snapshots
	DataSpace* DataspaceHVelY;   ///< dataspace for EACH Hole Velocity Y snapshots

	float *HDen;       // number density
	float *HVelX;      // fluid velocity x component
	float *HVelY;      // fluid velocity y component


	explicit DiracGraphene2D(SetUpParametersCNP &input_parameters);
		~DiracGraphene2D();

		/*!
		 * @brief Calculates @f$\Delta x@f$ and imposes Courant–Friedrichs–Lewy condition to @f$\Delta t@f$
		 *
		 * The method retrieves the spatial discretization @f$\Delta x = L / N_x @f$ and then sets the time step as @f$\Delta t = \Delta x/\lambda  @f$
		 * where @f$\lambda@f$ is given by
		 * @f{*eqnarray}
		   \lambda =1.2v_F \quad&\mathrm{if}\quad S<0.36v_F\\
		   \lambda =1.97S + 0.5 v_F \quad&\mathrm{otherwise} @f}
		 */
		void CflCondition() override;

		void InitialCondRand() override;             ///< Initial condition, zero velocity and constant density with 0.5% white noise
		virtual void InitialCondPointDen();					///< Initial condition, zero velocity and point of density with 0.5% white noise
        void InitialCondUniform();

		/*!
		 * @brief Sets the total simulation time in units of @f$v_0/L@f$
		 *
		 * The method predicts and set the maximum time for a typical Dyakonov-Shur instability simulation such that:
		  -# It is larger than the usual saturation time
		  -# After the saturation at least 10 oscillations occur

		  Our tests concluded heuristically that such criteria can be met setting @f$T_{max}=5+0.02S+20/S@f$
		 */
		void SetSimulationTime() override;

		/*!
		 * @brief Converts the mass density flux to velocity on the entire simulation grid.
		 *
		 * Since the mass of the fluid element is not a constant the in the graphene electronic fluid, one needs to perform the transformation
		   @f[ \vec{v} = \frac{\vec{p}}{n^{3/2}} @f]
		 * */
//		//void MassFluxToVelocity(string grid) override; // Converts the mass density flux back to velocity, in graphene  v = p n^{-3/2}
		/*Override fluxes and sources to specifics of graphene physics*/
		float EleDensitySource(StateVec2D Uelec , StateVec2D Uholes);   ///< density equation (continuity equation) source term
		float EleXMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< velocity X component equation (momentum equation) source term
		float EleYMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< velocity y component equation (momentum equation) source term

		float EleDensityFluxX(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< density equation (continuity equation) conserved flux X component
		float EleDensityFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes) ; ///< density equation (continuity equation) conserved1 flux Y component

		float EleXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes); ///< velocity X component equation (momentum equation) conserved flux X component
		float EleXMomentumFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< velocity X component equation (momentum equation) conserved flux Y component

		float EleYMomentumFluxX(StateVec2D Uelec , __attribute__((unused))  StateVec2D Uholes); ///< velocity Y component equation (momentum equation) conserved flux X component
		float EleYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes) ; ///< velocity Y component equation (momentum equation) conserved flux Y component

		float HolDensitySource(StateVec2D Uelec , StateVec2D Uholes);   ///< density equation (continuity equation) source term
		float HolXMomentumSource(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< velocity X component equation (momentum equation) source term
		float HolYMomentumSource(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes); ///< velocity y component equation (momentum equation) source term

		float HolDensityFluxX(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes); ///< density equation (continuity equation) conserved flux X component
		float HolDensityFluxY(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes); ///< density equation (continuity equation) conserved1 flux Y component

		float HolXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes); ///< velocity X component equation (momentum equation) conserved flux X component
		float HolXMomentumFluxY(__attribute__((unused))  StateVec2D Uelec , StateVec2D Uholes); ///< velocity X component equation (momentum equation) conserved flux Y component

		float HolYMomentumFluxX(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes); ///< velocity Y component equation (momentum equation) conserved flux X component
		float HolYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes); ///< velocity Y component equation (momentum equation) conserved flux Y component


	void RichtmyerStep1() override;
	void RichtmyerStep2() override;
		/*!
		 * @brief Writes the line of time t on the simplified .dat file output
		 *
		 * As a way to easily access the main results of the simulation the simplified .dat file stores the following quantities by columns:
		 * -# Time @f$t@f$
		 * -# Density at drain contact for electrons @f$n(x=L)@f$
		 * -# Mass flux along x at drain contact for electrons @f$p_x(x=L)@f$
		 * -# Density at source contact for electrons @f$n_x(x=0) @f$
		 * -# Mass flux along x at source contact for electrons @f$p_x(x=0) @f$
		 * -# Density at drain contact for holes @f$n(x=L)@f$
		 * -# Mass flux along x at drain contact for holes  @f$p_x(x=L)@f$
		 * -# Density at source contact for holes @f$n_x(x=0) @f$
		 * -# Mass flux along x at source contact for holes @f$p_x(x=0) @f$
		 * */
		void WriteFluidFile(float t) override; // writes the line of time t on the simplified .dat file output
		
		void VelocityLaplacianWeighted19() override;
	    void ForwardTimeOperator() override; ///< Time evolution for the FTCS method employed for the parabolic operators.
		void ForwardTimeOperator(char field) override; ///< Time evolution for the FTCS method employed for the parabolic operators.

	void CopyFields() override;
	void SaveSnapShot() override;
	void CreateHdf5File() override;

    std::ofstream phi_preview;
};


#endif //DIRACGRAPHENE2DLIBNOT_H