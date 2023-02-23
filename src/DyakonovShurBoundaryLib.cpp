/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/BoundaryLib.h"
#include "includes/DyakonovShurBoundaryLib.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// ONE-DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid1D& fluid_class) {
	int nx=fluid_class.SizeX();
	fluid_class.Umain[0].n()=1.0f;
	fluid_class.Umain[0].v()=fluid_class.Umain[1].v();
	fluid_class.Umain[nx-1].n()=fluid_class.Umain[nx-2].n();
	fluid_class.Umain[nx-1].v()=1.0f/fluid_class.Umain[nx-1].n();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// TWO-DIMENSIONAL FLUIDS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void DyakonovShurBoundaryCondition::DyakonovShurBc(GrapheneFluid2D& fluid_class) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();

	for(int j=0; j < ny; j++) {
/*
		fluid_class.Den[0 + j * nx] = 1.0f;            //constant density at x=0
		fluid_class.FlxX[0 + j * nx] =
				fluid_class.FlxX[1 + j * nx] * pow(fluid_class.Den[1 + j * nx], -1.5f);            //free flux at x=0
		fluid_class.FlxY[0 + j * nx] = 0.0f;                    //flux only on x at x=0

		fluid_class.Den[nx - 1 + j * nx] = fluid_class.Den[nx - 2 + j * nx];            //free density at x=L
		fluid_class.FlxX[nx - 1 + j * nx] = sqrt(fluid_class.Den[nx - 1 + j * nx]);    //constant current at x=L (flux equals mass)
		fluid_class.FlxY[nx - 1 + j * nx] = 0.0f;                    //idem at x=L
*/
		fluid_class.Umain[0 + j * nx].n() = 1.0f;            //constant density at x=0
		fluid_class.Umain[0 + j * nx].px() =
				fluid_class.Umain[1 + j * nx].px() * pow(fluid_class.Umain[1 + j * nx].n(), -1.5f);            //free flux at x=0
		fluid_class.Umain[0 + j * nx].py() = 0.0f;                    //flux only on x at x=0

		fluid_class.Umain[nx - 1 + j * nx].n() = fluid_class.Umain[nx - 2 + j * nx].n();            //free density at x=L
		fluid_class.Umain[nx - 1 + j * nx].px() = sqrt(fluid_class.Umain[nx - 1 + j * nx].n());    //constant current at x=L (flux equals mass)
		fluid_class.Umain[nx - 1 + j * nx].py() = 0.0f;                    //idem at x=L


	}
}

void DyakonovShurBoundaryCondition::DSFeedbackBc(GrapheneFluid2D &fluid_class, float gain) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();

	for(int j=0; j < ny; j++) {
		/*
		fluid_class.Den[0 + j * nx] = 1.0f;            //constant density at x=0
		fluid_class.FlxX[0 + j * nx] =
				fluid_class.FlxX[1 + j * nx] * pow(fluid_class.Den[1 + j * nx], -1.5f);            //free flux at x=0
		fluid_class.FlxY[0 + j * nx] = 0.0f;                    //flux only on x at x=0

		fluid_class.Den[nx - 1 + j * nx] = fluid_class.Den[nx - 2 + j * nx];            //free density at x=L
		fluid_class.FlxX[nx - 1 + j * nx] = sqrt(fluid_class.Den[nx - 1 + j * nx])
				+ gain*fluid_class.FlxX[0 + j * nx] ;    //constant current at x=L (flux equals mass) plus gain* current at x=0
		fluid_class.FlxY[nx - 1 + j * nx] = 0.0f;                    //idem at x=L
		*/


		fluid_class.Umain[0 + j * nx].n() = 1.0f;            //constant density at x=0
		fluid_class.Umain[0 + j * nx].px() =
				fluid_class.Umain[1 + j * nx].px() * pow(fluid_class.Umain[1 + j * nx].n(), -1.5f);            //free flux at x=0
		fluid_class.Umain[0 + j * nx].py() = 0.0f;                    //flux only on x at x=0

		fluid_class.Umain[nx - 1 + j * nx].n() = fluid_class.Umain[nx - 2 + j * nx].n();            //free density at x=L
		fluid_class.Umain[nx - 1 + j * nx].px() = sqrt(fluid_class.Umain[nx - 1 + j * nx].n())  + gain*fluid_class.Umain[0 + j * nx].px() ;    //constant current at x=L (flux equals mass) plus gain* current at x=0;
		fluid_class.Umain[nx - 1 + j * nx].py() = 0.0f;                    //idem at x=L

	}
}
