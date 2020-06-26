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

#include "Tethys1DLib.h"


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

using namespace H5;
using namespace std;
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);




int main(int argc, char **argv){

	const int Npoints=201; 							// number of spatial points
	float t=0.0;
	float dx;								// spatial discretisation
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
	
	
	GrapheneFluid1D	graph(Npoints,input_vel_snd, input_vel_fer, input_kin_vis,input_col_freq);
	

	
	graph.BannerDisplay();
	/*......CFL routine to determine dt...............................*/	
	graph.CFLCondition();
	dx=graph.GetDx();
	dt=graph.GetDt();
//	graph.SetTmax(10.0);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	graph.SetSimulationTime();
	float T_max=graph.GetTmax();
	/*................................................................*/

	/*.........Output files and streams...............................*/
	ElectroAnalysis elec;
	elec.CreateElectroFile(graph);
	graph.CreateFluidFile();
	graph.CreateHDF5File();
	
	/*................................................................*/
	
	graph.WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),graph.GetKinVis(), dt, dx, T_max);
	RecordLogFile(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, T_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	graph.BoundaryCond(3);
	////////////////////////////////////////////////////////////////////
	int time_step=0;
	int snapshot_per_Period = 10;   
	int points_per_Period = (2.0*MAT_PI/RealFreq(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),1))/dt;
	int snapshot_step = points_per_Period/snapshot_per_Period; 

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	//Main cycle
	while(t<=T_max && isfinite(graph.vel[(graph.SizeX()-1)/2])) // throw exception para nan / inf 
	{	
		++time_step;
		t += dt;
		// Main algorithm		
		graph.Richtmyer();
		// Impose boundary conditions
		graph.BoundaryCond(3);
		// Applying average filters for smoothing 	
		graph.Smooth(2);
		//Record full data
		if(data_save_mode && time_step % snapshot_step == 0 ){
			string str_time = to_string(time_step/snapshot_step);
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = graph.grp_den->createDataSet( name_dataset , hdf5_float, *graph.dataspace_den );
			dataset_den.write( graph.den_cor, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_vel = graph.grp_velX->createDataSet( name_dataset , hdf5_float, *graph.dataspace_velX );
			dataset_vel.write( graph.vel_cor, hdf5_float );
			dataset_vel.close();	
			
		}
		graph.WriteFluidFile(t);
		elec.WriteElectroFile(t,graph);
	}
	
	graph.WriteAtributes();
	
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	return 0;
}




