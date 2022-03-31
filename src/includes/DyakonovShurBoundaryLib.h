/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

/*!@file
 * @brief Header file for Dyakonov-Shur BC
 */

#ifndef DYAKONOVSHURBOUNDARYLIB_H
#define DYAKONOVSHURBOUNDARYLIB_H
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "includes/GrapheneFluid1DLib.h"


/*!
 * @brief Class for the asymmetric Dyakonov-Shur boundary conditions
 *
 * */
class  DyakonovShurBoundaryCondition : public DirichletBoundaryCondition {
public:
	/*!
	 * @brief Dyakonov-Shur boundary conditions 1D
	 *
	 * Implementation of the Dyakonov-Shur boundary conditions (1D case) of constant current at the drain and constant voltage (number density) at the source
	 * @f[ n(x=0)=n_0  @f] and @f[j(x=L)\equiv n(x=L)v(x=L)=n_0v_0 @f]
	 * */
	static void DyakonovShurBc(GrapheneFluid1D &fluid_class);

	/*!
	 * @brief Dyakonov-Shur boundary conditions 2D
	 *
	 * Implementation of the Dyakonov-Shur boundary conditions (2D case) of constant current at the drain and constant voltage (number density) at the source
	 * @f[ n(x=0)=n_0  @f] and @f[j(x=L)\equiv n(x=L)v(x=L)=n_0v_0 @f]
	 * */
	static void DyakonovShurBc(GrapheneFluid2D &fluid_class);

	/*!
	 * @brief Feedback Dyakonov-Shur boundary conditions 2D
	 *
	 * Implementation of the Dyakonov-Shur boundary conditions (2D case) with positive feedback at the drain side
	 * @f[ n(x=0)=n_0  @f] and @f[j(x=L)\equiv n(x=L)v(x=L)=n_0v_0 + \epsilon \times n(x=0)v(x=0) @f]
	 *
	 * @param gain feedback gain @f$\epsilon@f$
	 *
	 * */
	static void DSFeedbackBc(GrapheneFluid2D &fluid_class, float gain);
};


#endif //DYAKONOVSHURBOUNDARYLIB_H