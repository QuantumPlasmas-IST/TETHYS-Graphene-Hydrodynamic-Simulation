/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 2D Dirac-Fermi fluid
 */

#ifndef DIRACGRAPHENE2DLIB_H
#define DIRACGRAPHENE2DLIB_H

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

	float *hvel_snd_arr;    // array for saving the (potentially varying) S(x,y) function at main grid
	float *hvel_snd_arr_mid;// array for saving the (potentially varying) S(x,y) function at auxiliary grid
	float *hden_mid;       // mid or auxiliary grids defined with (Nx-1)*(Ny-1) size
	float *hflxX_mid;
	float *hflxY_mid;

	void ChooseGridPointers(const string &grid) override;

	float *hptr_den;
	float *hptr_px;
	float *hptr_py;

public :

	float *HDen;       // number density
	float *HVelX;      // fluid velocity x component
	float *HVelY;      // fluid velocity y component
	float *HFlxX;      // mass density flux x component
	float *HFlxY;      // mass density flux y component
	float *HCurX;      // current density x component
	float *HCurY;      // current density y component

	explicit DiracGraphene2D(SetUpParameters &input_parameters);
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
		float DensitySource(float n, float flx_x, float flx_y, float mass, float s)override;   ///< density equation (continuity equation) source term
		float XMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity X component equation (momentum equation) source term
		float YMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity y component equation (momentum equation) source term

		float DensityFluxX(GridPoint2D p, char side ) override; ///< density equation (continuity equation) conserved flux X component
		float DensityFluxY(GridPoint2D p, char side ) override; ///< density equation (continuity equation) conserved1 flux Y component

		float XMomentumFluxX(GridPoint2D p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux X component
		float XMomentumFluxY(GridPoint2D p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux Y component

		float YMomentumFluxX(GridPoint2D p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux X component
		float YMomentumFluxY(GridPoint2D p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux Y component


		float HDensityFluxX(GridPoint2D p, char side ); ///< density equation (continuity equation) conserved flux X component
		float HDensityFluxY(GridPoint2D p, char side ); ///< density equation (continuity equation) conserved1 flux Y component

		float HXMomentumFluxX(GridPoint2D p, char side ); ///< velocity X component equation (momentum equation) conserved flux X component
		float HXMomentumFluxY(GridPoint2D p, char side ); ///< velocity X component equation (momentum equation) conserved flux Y component

		float HYMomentumFluxX(GridPoint2D p, char side ); ///< velocity Y component equation (momentum equation) conserved flux X component
		float HYMomentumFluxY(GridPoint2D p, char side ); ///< velocity Y component equation (momentum equation) conserved flux Y component

		void Richtmyer() override;
};


#endif //DIRACGRAPHENE2DLIB_H