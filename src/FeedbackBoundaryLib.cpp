/************************************************************************************************\
* 2022 Pedro Cosme , João Santos , Ivan Figueiredo and Diogo Simões                              *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/BoundaryLib.h"
#include "includes/FeedbackBoundaryLib.h"

FeedbackBoundaryCondition::FeedbackBoundaryCondition(float time_delay, float delta){
    count=0;
    Nsteps=(int)(time_delay/delta);
    Curr = new float[Nsteps];
    Dens = new float[Nsteps];
    for(int i = 0; i < Nsteps; ++i){
        Curr[i]=0;
        Dens[i]=0;
    }
}

FeedbackBoundaryCondition::~FeedbackBoundaryCondition()=default;

void FeedbackBoundaryCondition::VoltageFeedbackBc(GrapheneFluid1D& fluid_class, float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    float Vi = fluid_class.Den[nx-1];
    float Ii = fluid_class.Vel[nx-1]*fluid_class.Den[nx-1];
    float Vf = Trans[0]*Vi+Trans[1]*Ii;
    float If = Trans[2]*Vi+Trans[3]*Ii;
    fluid_class.Den[0]+=Vf+intens*cos(omega*t);
    fluid_class.Vel[0]+=If/fluid_class.Den[0];
}

void FeedbackBoundaryCondition::CurrentFeedbackBc(GrapheneFluid1D& fluid_class, float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    float Vi = fluid_class.Den[nx-1]-1;
    float Ii = fluid_class.Vel[nx-1]*fluid_class.Den[nx-1]-1;
    float Vf = Trans[0]*Vi+Trans[1]*Ii;
    float If = Trans[2]*Vi+Trans[3]*Ii;
    fluid_class.Den[0]+=Vf;
    fluid_class.Vel[0]+=If/fluid_class.Den[0]+intens*cos(omega*t);
}

void FeedbackBoundaryCondition::VoltageDelayFeedbackBc(GrapheneFluid1D& fluid_class, float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    Dens[count%Nsteps] = fluid_class.Den[nx-1];
    Curr[count%Nsteps] = fluid_class.Vel[nx-1]*fluid_class.Den[nx-1];
    count=(count+1)%Nsteps;
    float Vf = Trans[0]*Dens[count]+Trans[1]*Curr[count];
    float If = Trans[2]*Dens[count]+Trans[3]*Curr[count];
    fluid_class.Den[0]+=Vf+intens*cos(omega*t);
    fluid_class.Vel[0]+=If/fluid_class.Den[0];
}

void FeedbackBoundaryCondition::CurrentDelayFeedbackBc(GrapheneFluid1D& fluid_class, float* Trans, float intens, float omega, float t) {
	int nx=fluid_class.SizeX();
    Dens[count%Nsteps] = fluid_class.Den[nx-1]-1;
    Curr[count%Nsteps] = fluid_class.Vel[nx-1]*fluid_class.Den[nx-1]-1;
    count=(count+1)%Nsteps;
    float Vf = Trans[0]*Dens[count]+Trans[1]*Curr[count];
    float If = Trans[2]*Dens[count]+Trans[3]*Curr[count];
    fluid_class.Den[0]+=Vf;
    fluid_class.Vel[0]+=If/fluid_class.Den[0]+intens*cos(omega*t);
}