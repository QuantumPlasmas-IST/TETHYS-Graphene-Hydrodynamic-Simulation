/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/BoundaryLib.h"
#include "includes/FeedbackBoundaryLib.h"

void FeedbackBoundaryCondition::DensityFeedbackBc(GrapheneFluid1D& fluid_class) {
	int nx=fluid_class.SizeX();
    fluid_class.Den[nx-1]=fluid_class.Den[nx-2];
    fluid_class.Den[0]=fluid_class.Den[nx-1];
    fluid_class.Vel[0]=fluid_class.Vel[1];
    fluid_class.Vel[nx-1]=fluid_class.Vel[nx-2];
}