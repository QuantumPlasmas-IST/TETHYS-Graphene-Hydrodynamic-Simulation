/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/



#include "includes/InitialConditionLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "TethysBaseLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

float ff_top(float x){
	return 0.2*x;
}

float ff_bottom(float x){
	return (40+20*cos(x/20.) );
}

int main(int argc, char **argv){

	int NX = 400;
    int NY = 200;
	Geometry Geom(NX,NY);
	Geom.fronteira.D.set_Domain(ff_top,ff_bottom);
    Geom.fronteira.set_Edge();
	Geom.dominio.dom = Geom.fronteira.D.dom;

	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	

	GrapheneFluid2D graph(parameters);

	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	//graph.SetDt(graph.GetDt()*0.25f);
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	float sound=graph.GetVelSnd();
	float ly=graph.GetLengthY();
	//graph.GetLengthX();

	//std::function<float(float,float)> variationS = [=](float x,float y){ return sound+5.f* tanh(6.0f*cos(2.0f*MAT_PI*2.0f*x)); };
	//std::function<float(float,float)> variationS = [=](float x,float y){ return sound*(1+0.3f*x* abs(y-0.5*ly)); };
	//graph.SetSound(variationS);
	graph.SetSound();
	//graph.SetSimulationTime();
	//graph.SetTmax(3.f);

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
//	graph.InitialCondGeneral([](float x,float y) { return 1.0f+0.1f/cosh(10.0f*sqrt((x-.5f)*(x-.5f)+(y-.5f)*(y-.5f))); },[](float x,float y) { return 0.5f; },[](float x,float y) { return 0.0f; });
//	graph.InitialCondGeneral([](float x,float y) { return 1.0f+0.1f/cosh(10.0f*(x-.5f)); },[](float x,float y) { return 0.5f/cosh(10.0f*(x-.5f)); },[](float x,float y) { return 0.0f; });
//	graph.InitialCondGeneral([](float x,float y) { return 0.8; },[](float x,float y) { return 0.5f; },[](float x,float y) { return 0.0f; });

	//graph.InitialCondRand();

InitialCondition::Rand(graph);

//	graph.InitialCondTest();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
	//BoundaryCondition::SetSlope(0.0f);
	//BoundaryCondition::SetBottomEdge(graph);
	//BoundaryCondition::SetTopEdge(graph);
	/*................................................................*/

	DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
	DirichletBoundaryCondition::YClosedFreeSlip(graph);

//	BoundaryCondition::XPeriodic(graph);
//	BoundaryCondition::YPeriodic(graph);


	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	while (t <= graph.GetTmax() ){
		int percentage=100*GrapheneFluid2D::TimeStepCounter/(graph.GetTmax()/dt);
		cout << percentage<<"%\033[?25l"; //prints the percentage of simulation completed

		t += dt;
		GrapheneFluid2D::TimeStepCounter++;



		graph.Richtmyer(Geom);


		/*+++++++++++++++++++++++++++++++++++++*
		 * Change the boundary conditions here *
		 *+++++++++++++++++++++++++++++++++++++*/
		DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		DirichletBoundaryCondition::YClosedNoSlip(graph,Geom);
		//DirichletBoundaryCondition::YClosedFreeSlip(graph);

//		BoundaryCondition::XPeriodic(graph);
//		BoundaryCondition::YPeriodic(graph);


	//	if(graph.GetThermDiff()!=0.0){
	//		DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f);  // 300K corresponds to 0.22*Fermi temperature
	//	}
		if(graph.GetKinVis()!=0.0f ) {
			graph.ParabolicOperatorWeightedExplicit19('V');
			//*+++++++++++++++++++++++++++++++++++++*
			// * Change the boundary conditions here *
			// *+++++++++++++++++++++++++++++++++++++
			DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
			DirichletBoundaryCondition::YClosedNoSlip(graph);
			//DirichletBoundaryCondition::YClosedFreeSlip(graph);
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