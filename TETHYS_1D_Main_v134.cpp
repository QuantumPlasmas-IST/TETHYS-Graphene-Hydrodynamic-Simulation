#include "Tethys1DLib.h"
#include "BoundaryLib.h"
#include "ElectricLib.h"
#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;






int main(int argc, char **argv){
	const int npoints=101; 							// number of spatial points
	float t=0.0;
	//float dx;								// spatial discretisation
	float dt;								// time step



	SetUpInput parameters(argc, argv);

	GrapheneFluid1D	graph(npoints, parameters);
	DyakonovShurBoundaryCondition BC;
	


	graph.BannerDisplay();
	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	//dx=graph.GetDx();
	dt=graph.GetDt();
//	graph.SetTmax(10.0);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	graph.SetSimulationTime();
	float t_max=graph.GetTmax();
	/*................................................................*/

	/*.........Output files and streams...............................*/
	ElectroAnalysis elec;
	elec.CreateElectroFile(graph);
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	/*................................................................*/

	graph.WelcomeScreen();
	//Record_Log_File(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, 0.0, t_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	BC.DyakonovShurBc(graph);
	////////////////////////////////////////////////////////////////////

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	//Main cycle
	while(t <= t_max ) {

		t += dt;
		graph.TimeStepCounter++;
		// Main algorithm		
		graph.Richtmyer();
		// Impose boundary conditions
		BC.DyakonovShurBc(graph);
		// Applying average filters for smoothing 	
		graph.Smooth(2);
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




