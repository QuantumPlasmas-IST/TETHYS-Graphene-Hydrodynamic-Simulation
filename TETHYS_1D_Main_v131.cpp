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


using namespace H5;
using namespace std;
const H5std_string   FILE_NAME( "Richtmyer_DATA_SET.h5" );
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);




int main(int argc, char **argv){
	/* Display name and version  */
    BannerDisplay();


	const int Npoints=201; 							// number of spatial points
	float t=0.0;
	float T_max=10.0;
//	const float leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
//	float vel_snd;						    // Sound speed
//	float vel_fer;							// Fermi velocity
//	float col_freq; 								// mean free path in units of GFET length



GrapheneFluid1D	graph(Npoints);

 	
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
	dt=graph.GetDt();
	graph.SetTmax(10.0);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	/*................................................................*/


	/*.........Output files and streams...............................*/
	graph.CreateElectroFile();
	graph.CreateFluidFile();
	/*................................................................*/

	
	
	/*.............  HDF5 file .......................................*/  	
	// throw HDF5 exceptions ???
	/*
	 * Create a file.
	 */
	H5File* hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	
	/*
	 * Create the groups ("folders") in the file
	 */
	Group* grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	Group* grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));
	Group* grp_vel = new Group( hdf5file->createGroup( "/Data/Velocity" ));
	Group* grp_cur = new Group( hdf5file->createGroup( "/Data/Current" ));
	/*
	 * Create attributes 
	 */	
	hsize_t dim_atr[1] = { 1 };
	// Create the data space for the attribute.
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = grp_dat->createAttribute( "S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = grp_dat->createAttribute( "Fermi velocity", hdf5_float, atr_dataspace);
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
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_total_time.write( hdf5_float, &T_max);
	atr_num_space_points.write( hdf5_int, &Npoints);
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
	
	const int RANK = 1; // 1D simulation
	hsize_t     dim_snap[1];              // dataset dimensions
	dim_snap[0] = Npoints; 
	DataSpace dataspace_den( RANK, dim_snap );
	DataSpace dataspace_vel( RANK, dim_snap );
	DataSpace dataspace_cur( RANK, dim_snap );

	/*................................................................*/
	
	WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),graph.GetKinVis(), dt, dx, T_max);
	RecordLogFile(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, T_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	graph.BoundaryCond(3);
	////////////////////////////////////////////////////////////////////
	
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	int time_step=0;
	
	while(t<=T_max && isfinite(graph.vel[(graph.SizeX()-1)/2])) // throw exception para nan / inf 
	{	
		++time_step;
		t += dt;
		
		graph.Richtmyer();
		
		// Impose boundary conditions
		graph.BoundaryCond(3);
		
		// Applying average filters for smoothing 	
		graph.Smooth(2);
		
		
		if(data_save_mode && time_step % 35 == 0 ){
		//Record full data
			string str_time = to_string(time_step/35);
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = grp_den->createDataSet( name_dataset , hdf5_float, dataspace_den );
			dataset_den.write( graph.den_cor, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_vel = grp_vel->createDataSet( name_dataset , hdf5_float, dataspace_vel );
			dataset_vel.write( graph.vel_cor, hdf5_float );
			dataset_vel.close();	
			
			DataSet dataset_cur = grp_cur->createDataSet( name_dataset , hdf5_float, dataspace_cur );
			dataset_cur.write( graph.cur_cor, hdf5_float );
			dataset_cur.close();
		}
		graph.WriteFluidFile(t);
		graph.WriteElectroFile(t);
		
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	atr_num_time_steps.write(hdf5_int, &time_step);
	atr_num_time_steps.close();


	grp_dat->close(); 
	grp_den->close(); 
	grp_vel->close();
	hdf5file->close();
	
	return 0;
}




