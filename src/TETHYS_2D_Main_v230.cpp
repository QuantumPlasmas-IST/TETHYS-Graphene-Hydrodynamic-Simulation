#include "includes/Tethys2DLib.h"
#include "includes/BoundaryLib.h"
#include "includes/ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){
	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	GrapheneFluid2D graph(parameters);
	DyakonovShurBoundaryCondition boundary_condition;
	//RobinBoundaryCondition boundary_condition;
	//DirichletBoundaryCondition boundary_condition;

	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	graph.SetTmax(0.5f);
	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHdf5File();



	if(parameters.SaveMode){
		graph.SaveSound();
	}
	/*................................................................*/

	GrapheneFluid2D::BannerDisplay();
	graph.WelcomeScreen();

	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	////////////////////////////////////////////////////////////////////
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;

	DyakonovShurBoundaryCondition::SetSlope(0.0f);
	DyakonovShurBoundaryCondition::SetBottomEdge(graph);
	DyakonovShurBoundaryCondition::SetTopEdge(graph);

	while (t <= graph.GetTmax() ){

		t += dt;
		GrapheneFluid2D::TimeStepCounter++;

		graph.Richtmyer();
		DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		//boundary_condition.YFree(graph);
		//boundary_condition.YClosedNoSlip(graph);
		DyakonovShurBoundaryCondition::YClosedFreeSlip(graph);
		/*
		boundary_condition.YClosedNoSlip(graph);
		boundary_condition.DensityLeft(graph, 1.0f);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.XFreeRight(graph);
*/

		if(graph.GetKinVis()!=0.0f || graph.GetCycFreq()!=0.0f) {
			//graph.ParabolicOperatorFtcs();
			graph.ParabolicOperatorWeightedExplicit19();
			DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
			//boundary_condition.YFree(graph);
			//boundary_condition.YClosedNoSlip(graph);
			DyakonovShurBoundaryCondition::YClosedFreeSlip(graph);


			/*
			boundary_condition.YClosedNoSlip(graph);
			boundary_condition.DensityLeft(graph, 1.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.XFreeRight(graph);
			*/
		}


		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}
		graph.WriteFluidFile(t); //TODO we caould probably save less points of the the fluid file
	}
	//Record atributes on hdf5 file
	if(parameters.SaveMode) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();


	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}