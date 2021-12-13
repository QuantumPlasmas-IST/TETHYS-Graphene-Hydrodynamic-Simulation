/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/



#include "includes/Fluid2DLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"

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

	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	graph.SetDt(graph.GetDt()*0.25f);
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	//graph.SetTmax(5.0f);

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

	/*...............Initialization...................................*/
	graph.InitialCondWave();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
	BoundaryCondition::SetSlope(0.0f);
	BoundaryCondition::SetBottomEdge(graph);
	BoundaryCondition::SetTopEdge(graph);
	/*................................................................*/


	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	while (t <=  graph.GetTmax() ){
		int percentage=100*GrapheneFluid2D::TimeStepCounter/(graph.GetTmax()/dt);
		cout << percentage<<"%\033[?25l"; //prints the percentage of simulation completed

		t += dt;
		GrapheneFluid2D::TimeStepCounter++;

		graph.Richtmyer();

		/*+++++++++++++++++++++++++++++++++++++*
		 * Change the boundary conditions here *
		 *+++++++++++++++++++++++++++++++++++++*/
		//DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		//DirichletBoundaryCondition::YClosedNoSlip(graph);

		BoundaryCondition::XPeriodic(graph);
		BoundaryCondition::YPeriodic(graph);

		if(graph.GetThermDiff()!=0.0){
			DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f);  // 300K corresponds to 0.22*Fermi temperature
		}

		if(graph.GetKinVis()!=0.0f || graph.GetThermDiff()!=0.0f  ) {
			graph.ParabolicOperatorWeightedExplicit19();

			/*+++++++++++++++++++++++++++++++++++++*
			 * Change the boundary conditions here *
			 *+++++++++++++++++++++++++++++++++++++*/
		//	DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		//	DirichletBoundaryCondition::YClosedNoSlip(graph);


			BoundaryCondition::XPeriodic(graph);
			BoundaryCondition::YPeriodic(graph);


			if(graph.GetThermDiff()!=0.0){
				DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f); // 300K corresponds to 0.22*Fermi temperature
			}
		}

		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}
		if( !( GrapheneFluid2D::TimeStepCounter % 2) ){
			graph.WriteFluidFile(t);
		}
		cout <<"\033[1G\033[2K"; //clears percentage of completion
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