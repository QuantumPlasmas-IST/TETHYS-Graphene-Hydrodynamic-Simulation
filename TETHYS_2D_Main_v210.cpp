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

const FloatType      HDF5FLOAT(PredType::NATIVE_FLOAT);


int main(int argc, char **argv){
	float t_max;
	t_max = 6;
	int npoints_x = 101;
	int npoints_y = 101;
	int npoints = npoints_x * npoints_y;
	
	float t=0.0;
	float dx,dy;								// spatial discretisation
	float dt;								// time step

	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis;
	Parameter_Initalization(argc, argv, data_save_mode, input_vel_snd, input_vel_fer, input_col_freq, input_kin_vis);
	
	
	GrapheneFluid2D	graph(npoints_x, npoints_y, input_vel_snd, input_vel_fer, input_kin_vis, input_col_freq);
	graph.BannerDisplay();
	BoundaryCondition::DyakonovShur BC;
	//BoundaryCondition::Dirichlet BCD;
	//BoundaryCondition BC;
	
	/*......CFL routine to determine dt...............................*/
	graph.CFLCondition();
	dx=graph.GetDx();
	dy=graph.GetDy();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/

	graph.SetSound();
	//graph.SetSimulationTime();
	graph.SetTmax(6.0);
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

	
	while(t <= t_max && isfinite(graph.VelX[npoints / 2])) // throw exception para nan / inf
	//while(time_step<=10000 && isfinite(graph.velX[npoints/2])) // throw exception para nan / inf
	{	
		++time_step;
		t += dt;

		//graph.SourceFTCS();
		// Impose boundary conditions
		//BC.X(graph);
		//BC.YFree(graph);

		graph.ViscosityFTCS();
		BC.X(graph);
		BC.YFree(graph);


		graph.Richtmyer();
		// Impose boundary conditions
		BC.X(graph);
		BC.YFree(graph);

//		graph.SourceFTCS();
		// Impose boundary conditions
//		BC.X(graph);
		//BC.YFree(graph);


		//graph.MagneticSource();
		//graph.SourceFTCS();

		graph.ViscosityFTCS();
		BC.X(graph);
		BC.YFree(graph);


//		BCD.Density(graph,1.0f,2.0f,1.0f,1.0f);
//		BCD.MassFluxX(graph,-1.0f,1.0f,0.0f,0.0f);

//		BC.YFree(graph);
		// Applying average filters for smoothing 	
		//graph.Smooth(2);

		
		if(data_save_mode && time_step % snapshot_step  == 0 ){
		//Record full data
			graph.MassFluxToVelocity(); //if placed here the profiling percentage drops from 7.7% to <1%
			string str_time = to_string(time_step/snapshot_step );
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = graph.GrpDen->createDataSet(name_dataset , HDF5FLOAT, *graph.DataspaceDen );
			dataset_den.write(graph.Den, HDF5FLOAT );
			dataset_den.close();
			
			DataSet dataset_vel_x = graph.GrpVelX->createDataSet(name_dataset , HDF5FLOAT, *graph.DataspaceVelX );
			dataset_vel_x.write(graph.VelX, HDF5FLOAT );
			dataset_vel_x.close();

			DataSet dataset_vel_y = graph.GrpVelY->createDataSet(name_dataset , HDF5FLOAT, *graph.DataspaceVelY );
			dataset_vel_y.write(graph.VelY, HDF5FLOAT );
			dataset_vel_y.close();
		}
		graph.WriteFluidFile(t);
		
		
	
		
	}
	

	if(data_save_mode ) {
		graph.WriteAtributes();
	}
	graph.CloseHDF5File();
	if(!data_save_mode ) {
		system("rm hdf5_2D*");
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;


	
	return 0;
}




