/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for Robin BC
 */

#ifndef ROBINBOUNDARYLIB_H
#define ROBINBOUNDARYLIB_H
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DirichletBoundaryLib.h"

/*!
 * @brief Class for Robin type boundary conditions.
 *
 * This type of boundary conditions are characterized by a linear combination of the value and the normal derivative of the quantities at the boundary, i.e.
 * @f[ a U +  b \frac{\partial U}{\partial \xi}= c\quad {\rm at}\quad \partial \Omega\quad {\rm with}\quad a,b,c\in\mathbb{R}  @f]
 *
 * In the context of the fluid simulation it is useful for setting the _slip length_ conditions at the interfaces, where the fluid velocity is set as
 *
 * @f[ v_x = \pm\ell \frac{\partial v_x}{\partial y} \quad v_y=0@f]
 * */
class  RobinBoundaryCondition : public DirichletBoundaryCondition
{
public:
	static void SlipLength(Fluid2D& fluid_class,float slip_length);        ///< slip length boundary condition Vx = +-l dVx/dy and Vy=0 at y=0 and y=W
	static void SlipLengthTop(Fluid2D& fluid_class,float slip_length);     ///< slip length boundary condition Vx = -l dVx/dy and Vy=0 at y=W
	static void SlipLengthBottom(Fluid2D& fluid_class,float slip_length);  ///< slip length boundary condition Vx = l dVx/dy and Vy=0 at y=0
};

#endif //ROBINBOUNDARYLIB_H
