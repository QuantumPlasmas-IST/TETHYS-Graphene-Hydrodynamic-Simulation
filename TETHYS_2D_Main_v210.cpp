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
	//int npoints = npoints_x * npoints_y;
	
	float t=0.0;
	float dx,dy;	// spatial discretisation
	float dt;		// time step



	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis,input_cyc_freq;
	Parameter_Initialization(argc, argv, data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis,
	                         input_cyc_freq);

	GrapheneFluid2D	graph(npoints_x, npoints_y, input_vel_snd, input_vel_fer, input_kin_vis, input_col_freq, input_cyc_freq);

	BoundaryCondition::DyakonovShur boundary_condition;
	BoundaryCondition::Dirichlet boundary_condition_Dirichelet;

	
	/*......CFL routine to determine dt...............................*/
	graph.SetLengthX(1.0f);
	graph.SetLengthY(2.0f);
	graph.CFLCondition();
	dx=graph.GetDx();
	dy=graph.GetDy();
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
	graph.CreateHDF5File();
	/*................................................................*/

	//t_max=3.0f; //encurtar o tempo para testes

	graph.BannerDisplay();
	graph.WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), graph.GetKinVis(), dt, dx,dy, t_max);
	Record_Log_File(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, dy, t_max);
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	////////////////////////////////////////////////////////////////////
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	int time_step=0;
	int snapshot_per_period = 10;
	int points_per_period = static_cast<int>((2.0 * MAT_PI /Real_Freq(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), 1)) / dt);
	int snapshot_step = points_per_period / snapshot_per_period;


	while (t <= t_max ){
		++time_step;
		t += dt;
		graph.Richtmyer();
		boundary_condition.X(graph);
		//boundary_condition_Dirichelet.MassFluxX(graph,1.0f,1.0f,0.0f,0.0f);
		//boundary_condition_Dirichelet.MassFluxY(graph,0.0f,0.0f,0.0f,0.0f);
		boundary_condition.YFree(graph);
		//boundary_condition.YClosedNoSlip(graph);


		/*if(graph.GetCycFreq()!=0.0f){
			graph.MagneticSourceFTCS();
			boundary_condition.X(graph);
			//boundary_condition_Dirichelet.MassFluxX(graph,1.0f,1.0f,0.0f,0.0f);
			//boundary_condition_Dirichelet.MassFluxY(graph,0.0f,0.0f,0.0f,0.0f);
			boundary_condition.YFree(graph);
			//boundary_condition.YClosedNoSlip(graph);

		}*/
		if(graph.GetKinVis()!=0.0f) {
			graph.ViscosityFTCS();
			boundary_condition.X(graph);
			//boundary_condition_Dirichelet.MassFluxX(graph,1.0f,1.0f,0.0f,0.0f);
			//boundary_condition_Dirichelet.MassFluxY(graph,0.0f,0.0f,0.0f,0.0f);
			boundary_condition.YFree(graph);
			//boundary_condition.YClosedNoSlip(graph);
		}

		//Record full hdf5 data
		if (data_save_mode && time_step % snapshot_step == 0) {
			graph.SaveSnapShot(time_step,snapshot_step);
		}
		graph.WriteFluidFile(t);
	}
	//Record atributes on hdf5 file
	if(data_save_mode ) {
		graph.WriteAtributes();
	}
	graph.CloseHDF5File();
	if(!data_save_mode ) {
		//Remove the empty hdf5 file if unused
		system("rm hdf5_2D*");
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}