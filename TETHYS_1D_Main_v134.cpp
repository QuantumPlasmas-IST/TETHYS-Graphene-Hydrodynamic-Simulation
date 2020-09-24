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
	float dx;								// spatial discretisation
	float dt;								// time step

	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis,input_cyc_freq=0.0;
	Parameter_Initialization(argc, argv, data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis,
	                         input_cyc_freq);

	
	GrapheneFluid1D	graph(npoints, input_vel_snd, input_vel_fer, input_kin_vis, input_col_freq);
	DyakonovShurBoundaryCondition BC;
	


	graph.BannerDisplay();
	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	dx=graph.GetDx();
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

	graph.WelcomeScreen(graph.GetVelFer());
	Record_Log_File(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, 0.0, t_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	BC.DyakonovShurBc(graph);
	////////////////////////////////////////////////////////////////////
	int time_step=0;
	int snapshot_per_period = 10;
	int points_per_period = static_cast<int>((2 * MAT_PI /
			Real_Freq(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), 1)) / dt);
	int snapshot_step = points_per_period / snapshot_per_period;

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	//Main cycle
	while(t <= t_max ) {
		++time_step;
		t += dt;
		// Main algorithm		
		graph.Richtmyer();
		// Impose boundary conditions
		BC.DyakonovShurBc(graph);
		// Applying average filters for smoothing 	
		graph.Smooth(2);
		//Record full data
		if(data_save_mode && time_step % snapshot_step == 0 ){
			graph.SaveSnapShot(time_step,snapshot_step);
		}
		graph.WriteFluidFile(t);
		elec.WriteElectroFile(t,graph);
	}
	if(data_save_mode ) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();
	if(!data_save_mode ) {
		system("rm hdf5_1D*");
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	return 0;
}




