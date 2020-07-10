#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   
#include <cassert>

#include <H5Cpp.h>

#include "Tethys2DLib.h"
#include "BoundaryLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


using namespace H5;
using namespace std;

const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);




int main(int argc, char **argv){
    float T_max=6;
	int NpointsX = 101;
	int NpointsY = 101;
	int Npoints = NpointsX*NpointsY;
	
	float t=0.0;
	float dx,dy;								// spatial discretisation
	float dt;								// time step

 	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis;
	ParameterInitalization(argc,argv,data_save_mode,input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis);
	
	
	GrapheneFluid2D	graph(NpointsX,NpointsY,input_vel_snd, input_vel_fer, input_kin_vis,input_col_freq);
	graph.BannerDisplay();
    BoundaryCondition::DyakonovShur BC;
    //BoundaryCondition::Dirichlet BCD;
    //BoundaryCondition BC;
	
	/*......CFL routine to determine dt...............................*/
// TODO Check CFL in order to have different space discretizations
	graph.CFLCondition();
	dx=graph.GetDx();
	dy=graph.GetDy();
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
// TODO review the set sound and set max time processes
	graph.SetSound();
	//graph.SetSimulationTime();
	//float T_max=graph.GetTmax();
	graph.SetTmax(6.0);
	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHDF5File();
	/*................................................................*/
	
	graph.WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),graph.GetKinVis(), dt, dx, T_max);
	RecordLogFile(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx,dy, T_max);

		
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	////////////////////////////////////////////////////////////////////
	
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	int time_step=0;
	int snapshot_per_Period = 10;   
	int points_per_Period = static_cast<int>(
            (2.0 * MAT_PI / RealFreq(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), 1)) / dt);
	int snapshot_step = points_per_Period/snapshot_per_Period; 

	
	while(t<=T_max && isfinite(graph.velX[Npoints/2])) // throw exception para nan / inf
	//while(time_step<=10000 && isfinite(graph.velX[Npoints/2])) // throw exception para nan / inf
	{	
// TODO
//  1) try to implement strang splittind instead of godunov splitting

		++time_step;
		t += dt;

        graph.SourceFTCS();
        // Impose boundary conditions
        BC.X(graph);
        BC.YFree(graph);


		graph.Richtmyer();
        // Impose boundary conditions
        BC.X(graph);
        BC.YFree(graph);

		//graph.MagneticSource();
        graph.SourceFTCS();
        // Impose boundary conditions
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
			
			DataSet dataset_den = graph.grp_den->createDataSet( name_dataset , hdf5_float, *graph.dataspace_den );
			dataset_den.write( graph.den, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_velX = graph.grp_velX->createDataSet( name_dataset , hdf5_float, *graph.dataspace_velX );
			dataset_velX.write( graph.velX, hdf5_float );
			dataset_velX.close();

			DataSet dataset_velY = graph.grp_velY->createDataSet( name_dataset , hdf5_float, *graph.dataspace_velY );
			dataset_velY.write( graph.velY, hdf5_float );
			dataset_velY.close();	
		}
		graph.WriteFluidFile(t);
		
		
	
		
	}
	
	
	graph.WriteAtributes();
    graph.CloseHDF5File();
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;


	
	return 0;
}



