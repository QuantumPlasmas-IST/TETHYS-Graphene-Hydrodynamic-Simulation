/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for 1D Dirac-Fermi fluid
 */


#ifndef GRAPHENEFLUID1DLIB_H
#define GRAPHENEFLUID1DLIB_H

#include "includes/Fluid1DLib.h"

/*!
 * @brief Graphene electronic fluid class in one dimension.
 *
 * The GrapheneFluid1D class describes a  fluid governed by the usual continuity and Cauchy momemtum equations, where the mass of the fluid element is *not* constant.
 * It overrides class Fluid1D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class GrapheneFluid1D : public Fluid1D{
	public :
		explicit GrapheneFluid1D(SetUpParameters &input_parameters);
		~GrapheneFluid1D();
		/*Override CFL condition to the case of graphene equations */
		void CflCondition() override;
		/*Override fluxes and sources to specifics of graphene physics*/
		float DensityFlux(float n,float v,__attribute__((unused)) float s) override;
		float VelocityFlux(float n,float v,float dv,float s) override;
		float DensitySource(__attribute__((unused)) float n,__attribute__((unused)) float v, __attribute__((unused)) float s) override;
		float VelocitySource(__attribute__((unused)) float n,float v,__attribute__((unused)) float s) override;
};

#endif //GRAPHENEFLUID1DLIB_H
