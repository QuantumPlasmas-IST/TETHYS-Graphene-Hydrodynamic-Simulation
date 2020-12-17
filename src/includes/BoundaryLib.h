/************************************************************************************************\
* Copyright (c) 2020 Pedro Cosme and João Santos                                                 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "TethysMathLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

//TODO review the entire boundary class

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
		static void SetTopEdge(Fluid2D &fluid_class);
		static void SetBottomEdge(Fluid2D &fluid_class);
		static void SetSlope(float boundary_slope);
		static float GetSlope() ;
		static void XFree(Fluid1D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		static void XFree(Fluid2D& fluid_class);           ///< open boundaries at x=0 and x=L for all variables and zero tangent velocity
		static void XFreeLeft(Fluid2D& fluid_class);       ///< open boundaries at x=0 for all variables and zero tangent velocity Vy=0
		static void XFreeRight(Fluid2D& fluid_class);      ///< open boundaries at x=L for all variables and zero tangent velocity Vy=0
		static void XPeriodic(Fluid1D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables  and zero tangent velocity
		static void XPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(x=0)=u(x=L) for all variables and zero tangent velocity
		static void YFree(Fluid2D& fluid_class);           ///< open boundaries at y=0 and y=W for all variables
		static void YFreeTop(Fluid2D& fluid_class);        ///< open boundaries at y=W for all variables
		static void YFreeBottom(Fluid2D& fluid_class);     ///< open boundaries at y=0 for all variables
		static void YPeriodic(Fluid2D& fluid_class);       ///< periodic boundaries u(y=0)=u(y=W) for all variables
		static void YClosedFreeSlip(Fluid2D& fluid_class); ///< zero flux across y=0 and y=W and free tangent velocity Vx
		static void YClosedNoSlip(Fluid2D& fluid_class);   ///< zero flux across y=0 and y=W and zero tangent velocity Vx=0
};

/*!
 * @brief Class for constant Dirichelet boundary conditions
 *
 * The class implements boundary conditions where the variables of the system are set to a constant value along the boundary. For the fluid this can translate in a variety of situations and physical scenarios.
 *
 * */
class  DirichletBoundaryCondition : public BoundaryCondition
{
	public: 	
	static void Density(Fluid1D& fluid_class, float left, float right);                                    ///< Fixed density at boundary n(x=0)=left and n(x=L)=right
	static void Density(Fluid2D& fluid_class, float left, float right, float top, float bottom);           ///< Fixed density at boundary n(x=0)=left, n(x=L)=right, n(y=0)=bottom, n(y=W)=top
	static void VelocityX(Fluid1D& fluid_class, float left, float right);                                  ///< Fixed Velocity at boundary V(x=0)=left and V(x=L)=right
	static void MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom);         ///< Fixed mass density flux x component at boundary px(x=0)=left, px(x=L)=right, px(y=0)=bottom, px(y=W)=top
	static void MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom);         ///< Fixed mass density flux y component at boundary py(x=0)=left, py(x=L)=right, py(y=0)=bottom, py(y=W)=top
	static void Jet(Fluid2D& fluid_class, float left, float left_width, float right, float right_width);   ///< Jet configuration i.e. fixed flux x component at a portion of given with around the center of the edges x=0 and x=L. Useful to study turbulence onset
	static void DensityRight(Fluid2D& fluid_class, float right);       ///< Fixed density at boundary n(x=L)=right
	static void MassFluxXRight(Fluid2D& fluid_class, float right);     ///< Fixed mass density flux X component at boundary px(x=L)=right
	static void MassFluxYRight(Fluid2D& fluid_class, float right);     ///< Fixed mass density flux Y component at boundary py(x=L)=right
	static void DensityLeft(Fluid2D& fluid_class, float left);         ///< Fixed density at boundary n(x=0)=left
	static void MassFluxXLeft(Fluid2D& fluid_class, float left);       ///< Fixed mass density flux X component at boundary px(x=0)=left
	static void MassFluxYLeft(Fluid2D& fluid_class, float left);       ///< Fixed mass density flux Y component at boundary py(x=0)=left
	static void DensityTop(Fluid2D& fluid_class, float top);           ///< Fixed density at boundary n(y=W)=top
	static void MassFluxXTop(Fluid2D& fluid_class, float top);         ///< Fixed mass density flux X component at boundary px(y=W)=top
	static void MassFluxYTop(Fluid2D& fluid_class, float top);         ///< Fixed mass density flux Y component at boundary py(y=W)=top
	static void DensityBottom(Fluid2D& fluid_class, float bottom);     ///< Fixed density at boundary n(y=0)=bottom
	static void MassFluxXBottom(Fluid2D& fluid_class, float bottom);   ///< Fixed mass density flux x component at boundary px(y=0)=bottom
	static void MassFluxYBottom(Fluid2D& fluid_class, float bottom);   ///< Fixed mass density flux Y component at boundary py(y=0)=bottom
};

/*!
 * @brief Class for the asymmetric Dyakonov-Shur boundary conditions
 *
 * */
class  DyakonovShurBoundaryCondition : public DirichletBoundaryCondition
{
public:
	static void DyakonovShurBc(GrapheneFluid1D& fluid_class);  ///< Dyakonov-Shur boundary conditions 1D n(0)=1 n(L)V(L)=1
	static void DyakonovShurBc(GrapheneFluid2D& fluid_class);  ///< Dyakonov-Shur boundary conditions 2D n(x=0)=1 n(x=L)Vx(x=L)=1 Vy(x=0)=0 Vy(L=0)=0
};

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

#endif

