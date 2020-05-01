#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   

#include <H5Cpp.h>

#include "TethysLib.h"

using namespace H5;
using namespace std;
const H5std_string   FILE_NAME( "RichtmyerCLASSTEST_2D.h5" );
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);

int main()
{
	// Definig the class
	int NpointsX = 201;
	int NpointsY = 201;
	int Npoints = NpointsX*NpointsY;
	Fluid2D basic(NpointsX, NpointsY, 1);
	basic.CFLCondition();
	float input_vel_snd = 50.0;
	float dt = basic.GetDt();
	float dx = basic.GetDx();
	float dy = basic.GetDy();
	basic.SetVelSnd(input_vel_snd);
	basic.SetSound();
	const float Tmax = 0.5;
	int data_save_mode =1;

	// HDF5
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
	Attribute atr_dx = grp_dat->createAttribute( "Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = grp_dat->createAttribute( "Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = grp_dat->createAttribute( "Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_time_steps = grp_dat->createAttribute( "Number of time steps", hdf5_int, atr_dataspace);
	Attribute atr_num_space_points = grp_dat->createAttribute( "Number of spatial points", hdf5_int, atr_dataspace);
	// Write the attribute data. 
	atr_vel_snd.write( hdf5_float, &input_vel_snd);
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_total_time.write( hdf5_float, &Tmax);
	atr_num_space_points.write( hdf5_int, &Npoints);
	atr_vel_snd.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();

	const int RANK = 2; 		// 2D simulation
	hsize_t     dim_snap[2];   	// dataset dimensions
	dim_snap[0] = NpointsX;
	dim_snap[1] = NpointsY;
	
	DataSpace dataspace_den( RANK, dim_snap );
	DataSpace dataspace_velX( RANK, dim_snap );
	DataSpace dataspace_velY( RANK, dim_snap );

	basic.InitialCondRand();
	//basic.BoundaryCond(2);

	float t=0;
	int time_step = 0;

	cout << "Chegamos ao while\n";
	while(t<= Tmax)
	{
		if(data_save_mode && time_step % 35 == 0 ){
		//Record full data
			cout << "Chegamos à iterada : "<< time_step << "\n";
			for(int i=0; i<NpointsX; i+=50){
				for (int j=0; j<NpointsY; j+=50)
				{
					cout << "(" << i << "," << j << ")  " <<
					"den = " << basic.den[i][j] << " | " <<
					"velX = " << basic.velX[i][j] << " | " <<
					"velY = " << basic.velY[i][j] << endl;
				}
			}
			string str_time = to_string(time_step/35);
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = grp_den->createDataSet( name_dataset , hdf5_float, dataspace_den );
			dataset_den.write( basic.den, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_velX = grp_velX->createDataSet( name_dataset , hdf5_float, dataspace_velX );
			dataset_velX.write( basic.velX, hdf5_float );
			dataset_velX.close();

			DataSet dataset_velY = grp_velY->createDataSet( name_dataset , hdf5_float, dataspace_velY );
			dataset_velY.write( basic.velY, hdf5_float );
			dataset_velY.close();	
		}

		basic.Richtmyer();
		basic.BoundaryCond(2);
		//basic.Smooth(2);

		t+=dt;
		time_step++;
	}

	atr_num_time_steps.write(hdf5_int, &time_step);
	atr_num_time_steps.close();

	grp_dat->close(); 
	grp_den->close(); 
	grp_velX->close();
	grp_velY->close();

	return 0;
}