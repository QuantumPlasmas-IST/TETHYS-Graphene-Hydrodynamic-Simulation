/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"

class InitialCondition {
public:
	static void Rand(Fluid1D& fluid_class);
	static void Rand(Fluid2D& fluid_class);
	static void Rand(DiracGraphene2D& fluid_class);
	static void Test(Fluid2D& fluid_class);
	static void Wave(Fluid2D& fluid_class);
	static void InitialCondPointDen(DiracGraphene2D &fluid_class);
	static void InitialCondUniform(DiracGraphene2D &fluid_class);
	static void InitialCondUniform(DiracGraphene2D &fluid_class, Geometry *Geom);
	static void InitialCondGeneral(Fluid2D& fluid_class, function<float(float, float)> fden, function<float(float, float)> fvx, function<float(float, float)> fvy);

	static void InitialCondTest(Fluid1D& fluid_class);     ///< Initial condition for testing and debugging
	static void InitialCondGeneral(Fluid1D& fluid_class,function<float(float)> fden, function<float(float)> fvx);


};


#endif //INITIALCONDITION_H
