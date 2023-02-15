/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 2D Dirac-Fermi fluid
 */

#ifndef GRAPHENECNPFLUID2DLIBNOT_H
#define GRAPHENECNPFLUID2DLIBNOT_H

#include "includes/FluidCNP2DLibNOT.h"

/*!
 * @brief Graphene electronic fluid class in two dimensions.
 *
 * The GrapheneCNP2D class describes a  fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is *not* constant.
 * It overrides class FluidCNP2D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class GrapheneCNP2D : public FluidCNP2D{
protected:
	/*!
	 * @brief Converts the number density to efective mass
	 *
	 * Since the mass of the fluid element is not a constant the in the graphene electronic fluid, one needs to perform the transformation
	   @f[ m^\star = n^{3/2} @f]
	 * */
	float DensityToMass(float density) override;
	float vmax_x, vmax_y;
	float time_counter, save_T_e, save_T_h;
	float aux_var;
	
public :
	explicit GrapheneCNP2D(SetUpParametersCNP &input_parameters);
		~GrapheneCNP2D();

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
		void MassFluxToVelocity() override; // Converts the mass density flux back to velocity, in graphene  v = p n^{-3/2}
		/*Override fluxes and sources to specifics of graphene physics*/
		float DensitySource(float n, float flx_x, float flx_y, float mass, float s)override;   ///< density equation (continuity equation) source term
		float TemperatureSource(float n, float flx_x, float flx_y, float den_grad_x, float den_grad_y, float mass, float s) override;   ///< density equation (continuity equation) source term
		float XMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity X component equation (momentum equation) source term
		float YMomentumSource(float n, float flx_x, float flx_y, float mass, float s)override; ///< velocity y component equation (momentum equation) source term

		float DensityFluxX(GridPoint p, char side ) override; ///< density equation (continuity equation) conserved flux X component
		float DensityFluxY(GridPoint p, char side ) override; ///< density equation (continuity equation) conserved1 flux Y component

		float XMomentumFluxX(GridPoint p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux X component
		float XMomentumFluxY(GridPoint p, char side ) override; ///< velocity X component equation (momentum equation) conserved flux Y component

		float YMomentumFluxX(GridPoint p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux X component
		float YMomentumFluxY(GridPoint p, char side ) override; ///< velocity Y component equation (momentum equation) conserved flux Y component

		//TVD methods for charge-neutral BLG simulation
		void Reconstruction(const float* n, const float* n_h, const float* cur_x, const float* cur_y, const float* cur_x_h, const float* cur_y_h, const float* e, const float* e_h
			, const float* p, const float* p_h, const float* t, const float* t_h);
		float Get_slope_x(const float* arr, int m, int flag);
		float Get_slope_y(const float* arr, int m, int flag);
		float HLLFlux(float FL, float FR, float UL, float UR, float UL_star, float UR_star, float SL, float SR, float SM);
		void ComputeFluxX();
		void ComputeFluxY();
		float GetTmp(float p, float den, float t0);
		void MUSCL_RK3();

		//BLG fluxes
		float DensityFluxX2species(float cur);
		float DensityFluxY2species(float cur);
		float XCurrentFluxX2species(float den, float cur, float p);
		float XCurrentFluxY2species(float den, float cur1, float cur2);
		float YCurrentFluxX2species(float den, float cur1, float cur2);
		float YCurrentFluxY2species(float den, float cur, float p);
		float EnergyFluxX2species(float den, float cur, float e, float p);
		float EnergyFluxY2species(float den, float cur, float e, float p);

		//Self-consistent potential calculation
		void GetPhi2species(const float *n_e, const float *n_h);

		//Initial conditions for two fluid model
		void InitialCondRand2species();
};


#endif //GRAPHENECNPFLUID2DLIBNOT_H
