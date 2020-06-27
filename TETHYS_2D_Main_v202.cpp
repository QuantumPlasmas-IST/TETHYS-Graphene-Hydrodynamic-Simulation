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
	
    float T_max=1.5;
	int NpointsX = 201;
	int NpointsY = 201;
	int Npoints = NpointsX*NpointsY;
	
	float t=0.0;
	float dx,dy;								// spatial discretisation
	float dt;								// time step



 	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis;
	if(argc==6){
		input_vel_snd = atof(argv[1]);
		input_vel_fer = atof(argv[2]);
		input_col_freq = atof(argv[3]);
		input_kin_vis = atof(argv[4]);
		assert(atoi(argv[5])==0 || atoi(argv[5])==1);
		data_save_mode = atoi(argv[5]);	// full data or light save option
		}
	else{
		cout << "Define S value: "; // throw exceptions if the velocities or frequency are negative or if S<Vf
		cin >> input_vel_snd;
		cout << "Define vF value: ";
		cin >> input_vel_fer; 
		cout << "Define kinetic viscosity: ";
		cin >> input_kin_vis;
		cout << "Define collision frequency: ";
		cin >> input_col_freq;
		cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
		cin >> data_save_mode;
		}
	
	
	GrapheneFluid2D	graph(NpointsX,NpointsY,input_vel_snd, input_vel_fer, input_kin_vis,input_col_freq);

	graph.BannerDisplay();
	
	
    BoundaryCondition::DyakonovShur BC;

	
	/*......CFL routine to determine dt...............................*/	
	graph.CFLCondition();
	dx=graph.GetDx();
	dy=graph.GetDy();
	dt=graph.GetDt();
	graph.SetTmax(1.5);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	//graph.SetSimulationTime();
	//float T_max=graph.GetTmax();
	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHDF5File();
	/*................................................................*/
	
	graph.WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),graph.GetKinVis(), dt, dx, T_max);
	RecordLogFile(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, T_max);

		
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
//	graph.BoundaryCond(1);
	////////////////////////////////////////////////////////////////////
	
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	int time_step=0;
	int snapshot_per_Period = 10;   
	int points_per_Period = (2.0*MAT_PI/RealFreq(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),1))/dt;
	int snapshot_step = points_per_Period/snapshot_per_Period; 

	
	while(t<=T_max && isfinite(graph.velX[Npoints/2])) // throw exception para nan / inf 
	{	
		++time_step;
		t += dt;
		
		graph.Richtmyer();
		graph.MassFluxToVelocity();
		// Impose boundary conditions
//		graph.BoundaryCond(1);
		
BC.X(graph);
BC.YFree(graph);		
		
		// Applying average filters for smoothing 	
		//graph.Smooth(2);
		
		if(data_save_mode && time_step % snapshot_step  == 0 ){
		//Record full data
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
	
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;


	
	return 0;
}




