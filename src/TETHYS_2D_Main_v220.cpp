#include "Tethys2DLib.h"
#include "BoundaryLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){

	SetUpInput parameters(argc, argv);

	int npoints_x=100;
	int npoints_y=100;
	float size_x=1.0f;
	float size_y=1.0f;
	float aspect_ratio = parameters.AspectRatio;

	if(aspect_ratio>1.0f){
		size_x=1.0f*aspect_ratio;
		size_y=1.0f;
		npoints_y=101;
		npoints_x= static_cast<int>( (npoints_y-1)*aspect_ratio)+1;
	}
	if(aspect_ratio==1.0f){
		size_x=1.0f;
		size_y=1.0f;
		npoints_x=101;
		npoints_y=101;
	}
	if(aspect_ratio<1.0f){
		size_x=1.0f;
		size_y=1.0f/aspect_ratio;
		npoints_x=101;
		npoints_y= static_cast<int>( (npoints_x - 1) / aspect_ratio ) + 1;
	}


	float t=0.0;
	float dt;		// time step


	GrapheneFluid2D	graph(npoints_x, npoints_y, parameters);
	DyakonovShurBoundaryCondition boundary_condition;
	//RobinBoundaryCondition boundary_condition;
	//DirichletBoundaryCondition boundary_condition;

	/*......CFL routine to determine dt...............................*/
	graph.SetLengthX(size_x);
	graph.SetLengthY(size_y);
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

		if(graph.GetCycFreq()!=0.0f){
			graph.MagneticSourceFtcs();
			boundary_condition.DyakonovShurBc(graph);
			boundary_condition.YFree(graph);
		}

/*
		boundary_condition.YFree(graph);
		boundary_condition.DensityLeft(graph, 1.0f);
		boundary_condition.DensityRight(graph, 1.0f);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.MassFluxXRight(graph, 1.0f);
		*/
		//boundary_condition.MassFluxXRight(graph, 1.0f);
		//boundary_condition.DensityRight(graph, 1.0f);

		/*
		boundary_condition.YFree(graph); //para menter a densidade free em y=0
		boundary_condition.XFreeRight(graph);
		boundary_condition.MassFluxXBottom(graph, 0.0f);
		boundary_condition.MassFluxYBottom(graph, 0.0f);
		boundary_condition.DensityLeft(graph, 1.0f);
		boundary_condition.DensityTop(graph, 1.0f);
		boundary_condition.MassFluxXTop(graph, 1.0f);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.MassFluxYLeft(graph, 0.0f);
		*/
/*
		boundary_condition.XFreeRight(graph);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.MassFluxYLeft(graph, 0.0f);
		boundary_condition.SlipLength(graph,1.5f);
*/

		/*if(graph.GetCycFreq()!=0.0f){
			graph.MagneticSourceFtcs();
			boundary_condition.YFree(graph);
			boundary_condition.DensityLeft(graph, 1.0f);
			boundary_condition.DensityRight(graph, 1.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.MassFluxXRight(graph, 1.0f);
			//boundary_condition.MassFluxXRight(graph, 1.0f);
			//boundary_condition.DensityRight(graph, 1.0f);
		//}
		//if(graph.GetKinVis()!=0.0f) {
		//	graph.ViscosityFtcs();
		//	boundary_condition.DyakonovShurBc(graph);
		//	boundary_condition.YFree(graph);

			/*
			boundary_condition.YFree(graph); //para menter a densidade free em y=0
			boundary_condition.XFreeRight(graph);
			boundary_condition.MassFluxXBottom(graph, 0.0f);
			boundary_condition.MassFluxYBottom(graph, 0.0f);
			boundary_condition.DensityTop(graph, 1.0f);
			boundary_condition.MassFluxXTop(graph, 1.0f);
			boundary_condition.DensityLeft(graph, 1.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.MassFluxYLeft(graph, 0.0f);
			*/

	/*
			boundary_condition.XFreeRight(graph);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.MassFluxYLeft(graph, 0.0f);
			boundary_condition.SlipLength(graph,1.5f);
	*/
		//}

		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
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