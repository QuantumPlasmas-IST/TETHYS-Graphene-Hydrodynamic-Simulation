#include "includes/Fluid1DLib.h"
#include "includes/ElectricLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid1DLib.h"

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
	//DyakonovShurBoundaryCondition boundary_condition;

	GrapheneFluid1D::BannerDisplay();
	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/

	float sound = graph.GetVelSnd();
//	std::function<float(float)> variationS = [=](float x){ return sound+.5f* tanh(6.0f*cos(2.0f*MAT_PI*2.0f*x)); };
//	graph.SetSound(variationS);
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

	// falta ajustar parâmetros físicos!
	float As = 1.0f;
	float X0 = 0.5f*graph.GetLengthX();
	float D0 = 1.0f;
	float cs = 1.0f;
	float lambda = 20.0f*graph.GetLengthX();

	std::function<float(float)> fden = [=](float x){ return  D0 + As/(cosh(lambda*(x-X0))); };
	std::function<float(float)> fvx  = [=](float x){ return cs; };

	graph.InitialCondGeneral(fden, fvx);
	//graph.InitialCondRand();
	//DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
	//graph.InitialCondTest();

	/*................................................................*/

	std::cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;


	graph.SetTmax(10.0f);
	//Main cycle
	while(t <= graph.GetTmax()) {
		t += dt;
		GrapheneFluid1D::TimeStepCounter++;
		// Main algorithm
		graph.Richtmyer();
		//graph.RungeKuttaTVD();

		// Impose boundary conditions
		//DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		BoundaryCondition::XPeriodic(graph);

		//graph.BohmOperator(0.015);
		//BoundaryCondition::XPeriodic(graph);

		// Applying average filters for smoothing 	
		//graph.Smooth(2);
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

	std::cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	std::cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	return 0;
}




