#include "Tethys2DLib.h"
#include "BoundaryLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){

	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();
	//int npoints_x=200;
	//int npoints_y=200;
	//float size_x=1.0f;
	//float size_y=1.0f;
	//float aspect_ratio = parameters.AspectRatio;




	float t=0.0;
	float dt;		// time step


	GrapheneFluid2D graph(parameters);
	DyakonovShurBoundaryCondition boundary_condition;
	//RobinBoundaryCondition boundary_condition;
	//DirichletBoundaryCondition boundary_condition;

//	cout<<graph.GetLengthX()<<endl;
//	cout<<graph.GetLengthY()<<endl;

	/*......CFL routine to determine dt...............................*/
	//graph.SetLengthX(size_x);
	//graph.SetLengthY(size_y);
	graph.CflCondition();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	graph.SetTmax(7.0);
	/*................................................................*/

	/*.........Output files and streams...............................*/
	ElectroAnalysis elec;
	elec.CreateElectroFile(graph);
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	if(parameters.SaveMode){
		graph.SaveSound();
	}
	/*................................................................*/

	graph.BannerDisplay();
	graph.WelcomeScreen();

	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	////////////////////////////////////////////////////////////////////
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;



	while (t <= graph.GetTmax() ){
		t += dt;
		graph.TimeStepCounter++;

		graph.Richtmyer();
		boundary_condition.DyakonovShurBc(graph);
		boundary_condition.YFree(graph);
		/*
		boundary_condition.YClosedNoSlip(graph);
		boundary_condition.DensityLeft(graph, 1.0f);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.XFreeRight(graph);*/

		/*if(graph.GetCycFreq()!=0.0f){
			graph.MagneticSourceFtcs();
			boundary_condition.YClosedNoSlip(graph);
			boundary_condition.DensityLeft(graph, 1.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.XFreeRight(graph);
		}*/
		if(graph.GetKinVis()!=0.0f) {
			graph.ViscosityFtcs();
			boundary_condition.DyakonovShurBc(graph);
			boundary_condition.YFree(graph);
		/*
			boundary_condition.YClosedNoSlip(graph);
			boundary_condition.DensityLeft(graph, 1.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.XFreeRight(graph);*/
		}
		

		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
			elec.WriteElectroFile(t,graph);
		}
		graph.WriteFluidFile(t);

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