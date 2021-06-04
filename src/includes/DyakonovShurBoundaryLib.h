//
// Created by pcosme on 23/12/2020.
//

#ifndef DYAKONOVSHURBOUNDARYLIB_H
#define DYAKONOVSHURBOUNDARYLIB_H
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiricheletBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "includes/GrapheneFluid1DLib.h"


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


#endif //DYAKONOVSHURBOUNDARYLIB_H
