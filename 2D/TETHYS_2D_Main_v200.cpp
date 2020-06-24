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


using namespace H5;
using namespace std;
const H5std_string   FILE_NAME( "Richtmyer_DATA_2D_SET.h5" );
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);




int main(int argc, char **argv){
	/* Display name and version  */
    BannerDisplay();
	int NpointsX = 201;
	int NpointsY = 201;
	int Npoints = NpointsX*NpointsY;
	GrapheneFluid2D graph(NpointsX, NpointsY, 1);

	float t=0.0;
	float T_max=10.0;
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
		graph.SetVelSnd(input_vel_snd);
		graph.SetVelFer(input_vel_fer);
		graph.SetKinVis(input_kin_vis);
		graph.SetColFreq(input_col_freq);
		}
	else{
		cout << "Define S value: "; // throw exceptions if the velocities or frequency are negative or if S<Vf
		cin >> input_vel_snd;
		graph.SetVelSnd(input_vel_snd);
		cout << "Define vF value: ";
		cin >> input_vel_fer; 
		graph.SetVelFer(input_vel_fer);
		cout << "Define kinetic viscosity: ";
		cin >> input_kin_vis;
		graph.SetColFreq(input_kin_vis);
		cout << "Define collision frequency: ";
		cin >> input_col_freq;
		graph.SetColFreq(input_col_freq);
		cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
		cin >> data_save_mode;
		}
	
	/*......CFL routine to determine dt...............................*/	
	graph.CFLCondition();
	dx=graph.GetDx();
	dy=graph.GetDy();
	dt=graph.GetDt();
	graph.SetTmax(10.0);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	/*................................................................*/
	
	/*.............  HDF5 file .......................................*/  	
	/*
	 * Create a file.
	 */
	H5File* hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	
	/*
	 * Create the groups ("folders") in the file
	 */
	Group* grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	Group* grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));
	Group* grp_velX = new Group( hdf5file->createGroup( "/Data/VelocityX" ));
	Group* grp_velY = new Group( hdf5file->createGroup( "/Data/VelocityY" ));
	/*
	 * Create attributes 
	 */	
	hsize_t dim_atr[1] = { 1 };
	// Create the data space for the attribute.
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = grp_dat->createAttribute( "S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = grp_dat->createAttribute( "Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_kin_vis = grp_dat->createAttribute( "Kinetic viscosity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = grp_dat->createAttribute( "Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = grp_dat->createAttribute( "Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = grp_dat->createAttribute( "Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = grp_dat->createAttribute( "Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_time_steps = grp_dat->createAttribute( "Number of time steps", hdf5_int, atr_dataspace);
	Attribute atr_num_space_points = grp_dat->createAttribute( "Number of spatial points", hdf5_int, atr_dataspace);
	// Write the attribute data. 
	atr_col_freq.write(hdf5_float, &input_col_freq);
	atr_vel_fer.write( hdf5_float, &input_vel_fer);
	atr_vel_snd.write( hdf5_float, &input_vel_snd);
	atr_kin_vis.write(hdf5_float, &input_kin_vis); 
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_total_time.write( hdf5_float, &T_max);
	atr_num_space_points.write( hdf5_int, &Npoints);
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
	
	const int RANK = 2; // 2D simulation
	hsize_t     dim_snap[2];              // dataset dimensions
	dim_snap[0] = NpointsX;
	dim_snap[1] = NpointsY;

	DataSpace dataspace_den( RANK, dim_snap );
	DataSpace dataspace_velX( RANK, dim_snap );
	DataSpace dataspace_velY( RANK, dim_snap );

	/*................................................................*/
	
	WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, T_max);
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
		graph.BoundaryCond(1);
		

		
		
		// Applying average filters for smoothing 	
		//graph.Smooth(2);
		
		if(data_save_mode && time_step % snapshot_step  == 0 ){
		//Record full data
			string str_time = to_string(time_step/snapshot_step );
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = grp_den->createDataSet( name_dataset , hdf5_float, dataspace_den );
			dataset_den.write( graph.den, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_velX = grp_velX->createDataSet( name_dataset , hdf5_float, dataspace_velX );
			dataset_velX.write( graph.velX, hdf5_float );
			dataset_velX.close();

			DataSet dataset_velY = grp_velY->createDataSet( name_dataset , hdf5_float, dataspace_velY );
			dataset_velY.write( graph.velY, hdf5_float );
			dataset_velY.close();	
		}
		graph.WriteFluidFile(t);
		
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	atr_num_time_steps.write(hdf5_int, &time_step);
	atr_num_time_steps.close();


	grp_dat->close(); 
	grp_den->close(); 
	grp_velX->close();
	grp_velY->close();
	hdf5file->close();
	
	return 0;
}




