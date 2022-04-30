/************************************************************************************************\
* 2022 Pedro Cosme , João Santos , Ivan Figueiredo and Diogo Simões                              *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

/*!@file
 * @brief Header file for Feedback BC
 */

#ifndef FEEDBACKBOUNDARYLIB_H
#define FEEDBACKBOUNDARYLIB_H
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "includes/GrapheneFluid1DLib.h"


/*!
 * @brief Class for the feedback effects boundary conditions
 *
 * */

class  FeedbackBoundaryCondition : public DirichletBoundaryCondition {
public:
	/*!
	 * @brief Feedback boundary conditions 1D
	 *
	 * Implementation of single transistor feedback boundary conditions (1D case) of free current and same voltage at source and drain
	 * @f[ n(x=0)=n(x=L)  @f]
	 * */
	static void VoltageFeedbackBc(GrapheneFluid1D &fluid_class, float* Trans, float intens, float omega, float t);
    static void CurrentFeedbackBc(GrapheneFluid1D &fluid_class, float* Trans, float intens, float omega, float t);

    FeedbackBoundaryCondition(float time_delay, float delta);
    ~FeedbackBoundaryCondition();
    void VoltageDelayFeedbackBc(GrapheneFluid1D &fluid_class, float* Trans, float intens, float omega, float t);
    void CurrentDelayFeedbackBc(GrapheneFluid1D &fluid_class, float* Trans, float intens, float omega, float t);

protected:
    int count;
    int Nsteps;
    float* Curr;
    float* Dens;
};


#endif //FEEDBACKBOUNDARYLIB_H