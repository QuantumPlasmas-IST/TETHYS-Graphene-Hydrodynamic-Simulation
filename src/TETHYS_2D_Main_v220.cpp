#include "Tethys2DLib.h"
#include "BoundaryLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){
	float t_max;
	t_max = 6;
	int npoints_x = 101;
	int npoints_y = 201;
	float t=0.0;
	float dt;		// time step

	SetUpInput parameters(argc, argv);
	GrapheneFluid2D	graph(npoints_x, npoints_y, parameters);
	DyakonovShurBoundaryCondition boundary_condition;
	//RobinBoundaryCondition boundary_condition;
	//DirichletBoundaryCondition boundary_condition;

	/*......CFL routine to determine dt...............................*/
	graph.SetLengthX(1.0f);
	graph.SetLengthY(2.0f);
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

	//t_max=3.0f; //encurtar o tempo para testes

	graph.BannerDisplay();
	graph.WelcomeScreen();
//	Record_Log_File(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, dy, t_max);
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	////////////////////////////////////////////////////////////////////
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;

//t_max=4.0f;

	while (t <= t_max ){
		t += dt;
		graph.TimeStepCounter++;

		graph.Richtmyer();
		boundary_condition.DyakonovShurBc(graph);
		boundary_condition.YFree(graph);
		/*
		boundary_condition.YFree(graph);
		boundary_condition.DensityLeft(graph, 2.0f);
		boundary_condition.MassFluxXLeft(graph, 1.0f);
		boundary_condition.MassFluxXRight(graph, 1.0f);
		boundary_condition.DensityRight(graph, 1.0f);
		*/
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
			boundary_condition.DensityLeft(graph, 2.0f);
			boundary_condition.MassFluxXLeft(graph, 1.0f);
			boundary_condition.MassFluxXRight(graph, 1.0f);
			boundary_condition.DensityRight(graph, 1.0f);
		}*/
		if(graph.GetKinVis()!=0.0f) {
			graph.ViscosityFtcs();
			boundary_condition.DyakonovShurBc(graph);
			boundary_condition.YFree(graph);

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
		}

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