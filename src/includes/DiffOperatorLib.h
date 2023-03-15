//
// Created by pcosme on 14/03/23.
//

#ifndef DIFFOPERATORLIB_H
#define DIFFOPERATORLIB_H

#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"


class DiffOperator{
private:
	static void VelocityXGradient_bulk(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityXGradient_top(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityXGradient_bottom(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityXGradient_left(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityXGradient_right(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityXGradient_corners(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);

	static void VelocityYGradient_bulk(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityYGradient_top(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityYGradient_bottom(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityYGradient_left(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityYGradient_right(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
	static void VelocityYGradient_corners(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);

//	Fluid2D *fluid_pointer;

public:

//	DiffOperator(Fluid2D &fluid);
//	~DiffOperator()=default;

	static void VelocityGradient(Fluid2D &fluid_class,StateVec2D *Uarray, int size_x, int size_y);
};


#endif //DIFFOPERATORLIB_H
