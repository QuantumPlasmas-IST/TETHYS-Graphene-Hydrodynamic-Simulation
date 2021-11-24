/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

/*!@file
 * @brief Header file for BC base class
 */

#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"


/*Base class for general boundary conditions*/

/*!
 * @brief Base class for boundary conditions
 *
 * The class implements the more general and simple boundary conditions such as free flow or periodic ones.
 *
 * */
class BoundaryCondition {
	protected:
		static int * BottomEdge;
		static int * TopEdge;
		static float  Slope;
	public :
		static void SetTopEdge(Fluid2D &fluid_class); ///< Defines which grid points belong to the upper edge, i.e. y=W, useful for non-rectangular domains
		static void SetBottomEdge(Fluid2D &fluid_class); ///<  Defines which grid points belong to the lower edge, i.e. y=0, useful for non-rectangular domains
		static void SetSlope(float boundary_slope); ///< Sets the slope of the lateral edges, for the scenarion of non rectangular domains
		static float GetSlope() ; ///< Returns the slope of the lateral edges, for the scenarion of non rectangular domains
		static void XFree(Fluid1D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		static void XFree(Fluid2D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		/*!
		* @brief open boundaries for all variables along designated x edge.
		* Implemnts open boundaries for all quantities, at x=0 (if x_limit=0) or x=L (if x_limit=0) exclusively.
		* @param x_limit label to choose the edge
		* */
		static void XFree(Fluid2D &fluid_class, int x_limit);       ///< open boundaries at x=0 for all variables and zero tangent velocity Vy=0
		static void XPeriodic(Fluid1D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables  and zero tangent velocity
		static void XPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables and zero tangent velocity
		static void YFree(Fluid2D& fluid_class);           ///< open boundaries at y=0 and y=W for all variables
		/*!
		 * @brief open boundaries for all variables along designated y edge.
		 * Implemnts open boundaries for all quantities, at y=0 (if y_limit=0) or y=W (if y_limit=0) exclusively.
		 * @param y_limit label to choose the edge
		 * */
		static void YFree(Fluid2D &fluid_class, int y_limit);        ///< open boundaries at y=W for all variables
		static void YPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(y=0)=u(y=W) for all variables
		static void YClosedFreeSlip(Fluid2D& fluid_class); ///< zero flux across y=0 and y=W and free tangent velocity Vx
		static void YClosedNoSlip(Fluid2D& fluid_class);   ///< zero flux across y=0 and y=W and zero tangent velocity Vx=0

		static void XFreeLeft(Fluid2D &fluid_class); ///< open boundaries at x=0 for all variables
		static void XFreeRight(Fluid2D &fluid_class); ///< open boundaries at x=L for all variables
		static void YFreeTop(Fluid2D &fluid_class); ///< open boundaries at y=W for all variables
		static void YFreeBottom(Fluid2D &fluid_class); ///< open boundaries at y=0 for all variables
};

#endif

