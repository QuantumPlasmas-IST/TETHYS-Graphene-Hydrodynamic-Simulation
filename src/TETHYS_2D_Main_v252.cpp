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
	//parameters.PrintParameters();
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	GrapheneFluid2D graph(parameters);

	//graph.SetThermDiff(0.6);


	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	//graph.SetTmax(7.0f);

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
	graph.InitialCondRand();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
	DyakonovShurBoundaryCondition::SetSlope(0.0f);
	DyakonovShurBoundaryCondition::SetBottomEdge(graph);
	DyakonovShurBoundaryCondition::SetTopEdge(graph);
	/*................................................................*/


	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	while (t <=  graph.GetTmax() ){
		int percentage=100*GrapheneFluid2D::TimeStepCounter/(graph.GetTmax()/dt);
		cout << percentage<<"%\033[?25l";

		t += dt;
		//float forcing=1.0f + 0.2f*sin(t*64.0f);


		GrapheneFluid2D::TimeStepCounter++;

		graph.Richtmyer();
		DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		//DyakonovShurBoundaryCondition::YFree(graph);
//		DyakonovShurBoundaryCondition::YClosedFreeSlip(graph);
		DirichletBoundaryCondition::YClosedNoSlip(graph);

		if(graph.GetThermDiff()!=0.0){
			DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f);
		}

		//BoundaryCondition::YFreeTop(graph);
		/*
		BoundaryCondition::XFreeRight(graph);
		DirichletBoundaryCondition::DensityLeft(graph, 1.0f);
		DirichletBoundaryCondition::MassFluxXLeft(graph, forcing);
		DirichletBoundaryCondition::MassFluxYLeft(graph, 0.0f);
		DirichletBoundaryCondition::MassFluxYBottom(graph, 0.0f);
		DirichletBoundaryCondition::MassFluxXBottom(graph, 0.0f);

		DirichletBoundaryCondition::MassFluxYTop(graph, 0.0f);
		DirichletBoundaryCondition::MassFluxXTop(graph, 0.0f);
*/
		//RobinBoundaryCondition::SlipLengthBottom(graph, 1.5f);

		if(graph.GetKinVis()!=0.0f || graph.GetThermDiff()!=0.0f  ) {
			graph.ParabolicOperatorWeightedExplicit19();
			DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
//			DyakonovShurBoundaryCondition::YClosedFreeSlip(graph);
			//DyakonovShurBoundaryCondition::YFree(graph);
			DirichletBoundaryCondition::YClosedNoSlip(graph);
			if(graph.GetThermDiff()!=0.0){
				DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f);
			}
            //DirichletBoundaryCondition::Temperature(graph, 0.22f, 0.22f, 0.22f, 0.22f); // 300Kelvin sao aproximadamente 0.2 Temperatura de Fermi

/*
            //BoundaryCondition::YFreeTop(graph);
            BoundaryCondition::XFreeRight(graph);


            DirichletBoundaryCondition::DensityLeft(graph, 1.0f);
            DirichletBoundaryCondition::MassFluxXLeft(graph, forcing);
            DirichletBoundaryCondition::MassFluxYLeft(graph, 0.0f);
            DirichletBoundaryCondition::MassFluxYBottom(graph, 0.0f);
            DirichletBoundaryCondition::MassFluxXBottom(graph, 0.0f);

			DirichletBoundaryCondition::MassFluxYTop(graph, 0.0f);
			DirichletBoundaryCondition::MassFluxXTop(graph, 0.0f);

			//RobinBoundaryCondition::SlipLengthBottom(graph, 1.5f);*/
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