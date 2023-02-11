#include "includes/Fluid1DLib.h"
#include "includes/ElectricLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid1DLib.h"
#include "includes/FeedbackBoundaryLib.h"

#include <functional>

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){
	//const int npoints=101; 							// number of spatial points
	float t=0.0;
	//float dx;								// spatial discretisation
	float dt;								// time step

	SetUpParameters parameters(argc, argv);
	GrapheneFluid1D graph(parameters);

	GrapheneFluid1D::BannerDisplay();
	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/

	float sound = graph.GetVelSnd();
	graph.SetSound();
	graph.SetSimulationTime();

	/*................................................................*/

	/*.........Output files and streams...............................*/
	ElectroAnalysis elec;
	elec.CreateElectroFile(graph);
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	/*................................................................*/

	graph.WelcomeScreen();



	/*...............Initialization...................................*/
	graph.InitialCondRand();
	//DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
	//graph.InitialCondTest();

	/*................................................................*/

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;

	//In the case of feedbacked boundary conditions
	/*
    float* Hfeed = new float[4];
    Hfeed[0]=0.2;
    Hfeed[1]=0.0;
    Hfeed[2]=0.0;
    Hfeed[3]=0.0;
    FeedbackBoundaryCondition feed(0.26,dt);
    */

	graph.SetTmax(10.0);
	
    //Main cycle
	while(t <= graph.GetTmax() ) {
		t += dt;
		GrapheneFluid1D::TimeStepCounter++;
		
        // Main algorithm		
		graph.Richtmyer();

		// Impose boundary conditions
        DyakonovShurBoundaryCondition::DyakonovShurBc(graph);

		//In the case of feedback aply:
		/*DirichletBoundaryCondition::DensityLeft(graph,1);
        //DirichletBoundaryCondition::VelocityXLeft(graph,1);
        //BoundaryCondition::XFreeRight(graph);
        //feed.VoltageDelayFeedbackBc(graph,Hfeed,0.1,5*2*M_PI,t);
        */

		//Record full data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}
		graph.WriteFluidFile(t);
		elec.WriteElectroFile(t,graph);
	}
	if(parameters.SaveMode ) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();

	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	return 0;
}




