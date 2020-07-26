#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>


#include <H5Cpp.h>

#include "Tethys2DLib.h"
#include "BoundaryLib.h"


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace H5;
using namespace std;

const FloatType HDF5FLOAT(PredType::NATIVE_FLOAT);


int main(int argc, char **argv){
	float t_max;
	t_max = 6;
	int npoints_x = 101;
	int npoints_y = 101;
	//int npoints = npoints_x * npoints_y;
	
	float t=0.0;
	float dx,dy;	// spatial discretisation
	float dt;		// time step

	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis,input_cyc_freq;
	Parameter_Initalization(argc, argv, data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis,input_cyc_freq);

	GrapheneFluid2D	graph(npoints_x, npoints_y, input_vel_snd, input_vel_fer, input_kin_vis, input_col_freq, input_cyc_freq);
	graph.BannerDisplay();
	BoundaryCondition::DyakonovShur boundary_condition;
	BoundaryCondition::Dirichlet boundary_condition_Dirichelet;
	//BoundaryCondition boundary_condition;
	
	/*......CFL routine to determine dt...............................*/
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
	graph.CreateFluidFile();
	graph.CreateHDF5File();
	/*................................................................*/
	
	graph.WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), graph.GetKinVis(), dt, dx, t_max);
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
		//boundary_condition.YFree(graph);
//		boundary_condition_Dirichelet.Density(graph,1.0f,1.0f,1.0f,1.0f);
		//boundary_condition_Dirichelet.MassFluxX(graph,-1.0f,-1.0f,0.0f,0.0f);
		//boundary_condition_Dirichelet.Jet(graph, -1.0f, 0.8f, -1.0f, 0.8f);
		//boundary_condition.YFree(graph);
		boundary_condition.YFree(graph);

//		graph.MagneticSourceFTCS();
//		boundary_condition.X(graph);
//		boundary_condition.YFree(graph);
//		boundary_condition_Dirichelet.Density(graph,1.0f,1.0f,1.0f,1.0f);
//      boundary_condition_Dirichelet.MassFluxX(graph,-1.0f,-1.0f,0.0f,0.0f);
//	    boundary_condition_Dirichelet.Jet(graph, 0.0, 0.0, -1.0f, 0.3);
//		boundary_condition.YFree(graph);
		if(graph.GetCycFreq()!=0.0f){
			graph.MagneticSourceFTCS();
			boundary_condition.X(graph);
			boundary_condition.YFree(graph);
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
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}