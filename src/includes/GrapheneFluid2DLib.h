/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 2D Dirac-Fermi fluid
 */

#ifndef GRAPHENEFLUID2DLIB_H
#define GRAPHENEFLUID2DLIB_H

#include "includes/Fluid2DLib.h"



/*!
 * @brief Graphene electronic fluid class in two dimensions.
 *
 * The GrapheneFluid2D class describes a  fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is *not* constant.
 * It overrides class Fluid2D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class GrapheneFluid2D : public Fluid2D{
protected:
	/*!
	 * @brief Converts the number density to efective mass
	 *
	 * Since the mass of the fluid element is not a constant the in the graphene electronic fluid, one needs to perform the transformation
	   @f[ m^\star = n^{3/2} @f]
	 * */
	float DensityToMass(float density) override;

	//float coefBohm = 0.0001f;

public :
	explicit GrapheneFluid2D(SetUpParameters &input_parameters);
		~GrapheneFluid2D();
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



	float DensitySource(StateVec2D U)override;   ///< density equation (continuity equation) source term
	float TemperatureSource(StateVec2D U) override;   ///< density equation (continuity equation) source term
	float XMomentumSource(StateVec2D U)override; ///< velocity X component equation (momentum equation) source term
	float YMomentumSource(StateVec2D U)override; ///< velocity y component equation (momentum equation) source term

	float DensityFluxX(StateVec2D U ) override; ///< density equation (continuity equation) conserved flux X component
	float DensityFluxY(StateVec2D U ) override; ///< density equation (continuity equation) conserved1 flux Y component

	float XMomentumFluxX(StateVec2D U ) override; ///< velocity X component equation (momentum equation) conserved flux X component
	float XMomentumFluxY(StateVec2D U ) override; ///< velocity X component equation (momentum equation) conserved flux Y component

	float YMomentumFluxX(StateVec2D U) override; ///< velocity Y component equation (momentum equation) conserved flux X component
	float YMomentumFluxY(StateVec2D U ) override; ///< velocity Y component equation (momentum equation) conserved flux Y component

};


#endif //GRAPHENEFLUID2DLIB_H
