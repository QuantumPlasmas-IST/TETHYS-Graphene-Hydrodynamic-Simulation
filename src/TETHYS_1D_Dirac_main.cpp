#include "includes/Fluid1DLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/DiracGraphene1DLib.h"

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

	SetUpParametersCNP parameters(argc, argv);
	DiracGraphene1D graph(parameters);

	DiracGraphene1D::BannerDisplay();
	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/

//	float sound = graph.GetVelSnd();
//	std::function<float(float)> variationS = [=](float x){ return sound+.5f* tanh(6.0f*cos(2.0f*MAT_PI*2.0f*x)); };
//	graph.SetSound(variationS);
	graph.SetSound();
	graph.SetSimulationTime();

	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	/*................................................................*/

	graph.WelcomeScreen();


	/*...............Initialization...................................*/
	//graph.InitialCondRand();
	graph.InitialCondTest();

	/*................................................................*/

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;


	graph.SetTmax(10.0);
	//Main cycle
	while(t <= graph.GetTmax() ) {
		t += dt;
		GrapheneFluid1D::TimeStepCounter++;
		// Main algorithm		
		graph.Richtmyer();

		// Impose boundary conditions
		BoundaryCondition::XPeriodic(graph);

		//Record full data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}
		graph.WriteFluidFile(t);
	}
	if(parameters.SaveMode ) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();

	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	return 0;
}




