/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  This example writes a dataset to a new HDF5 file.
 */




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


using namespace H5;
using namespace std;

const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);

const H5std_string  FILE_NAME( "SDS.h5" );

const int   NX = 201;                    // dataset dimensions
const int   NY = 201;
const int   RANK = 2;
int main (void)
{
	/*
	* Data initialization.
	*/
	int i, j, Time;
	float data[NX][NY];          // buffer for data to write

	/*
	* Create a new file using H5F_ACC_TRUNC access,
	* default file creation properties, and default file
	* access properties.
	*/
	//H5File file( FILE_NAME, H5F_ACC_TRUNC );
	H5File* hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

	/*
	* Create the groups ("folders") in the file
	*/
	Group* grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	Group* grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));
	Group* grp_vel = new Group( hdf5file->createGroup( "/Data/Velocity" ));
	Group* grp_cur = new Group( hdf5file->createGroup( "/Data/Current" ));

	/*
	* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	hsize_t     dimsf[2];              // dataset dimensions
	dimsf[0] = NX;
	dimsf[1] = NY;
	DataSpace dataspace_den( RANK, dimsf );
	DataSpace dataspace_vel( RANK, dimsf );
	DataSpace dataspace_cur( RANK, dimsf );

	for(Time=1;Time <= 1000; ++Time){
		cout << Time<<endl;
		string str_time = to_string(Time);
		string name_dataset = "snapshot_"+str_time;

		for (j = 0; j < NX; j++){
			for (i = 0; i < NY; i++){
			data[j][i] = sin(6.0*i/200.0-4.0*j/200.0-Time*0.01) ;
			}
		}
		/*
		* Create a new dataset within the file using defined dataspace and
		* datatype and default dataset creation properties.
		*/
		DataSet dataset_den = grp_den->createDataSet( name_dataset, hdf5_float, dataspace_den ); 
		DataSet dataset_vel = grp_vel->createDataSet( name_dataset, hdf5_float, dataspace_vel ); 
		DataSet dataset_cur = grp_cur->createDataSet( name_dataset, hdf5_float, dataspace_cur ); 
		/*
		* Write the data to the dataset using default memory space, file
		* space, and transfer properties.
		*/
		dataset_den.write( data, hdf5_float );
		dataset_den.close();
		dataset_vel.write( data, hdf5_float );
		dataset_vel.close();
		dataset_cur.write( data, hdf5_float );
		dataset_cur.close();
	}
	grp_dat->close(); 
	grp_den->close(); 
	grp_vel->close();
	hdf5file->close();
	return 0;  // successfully terminated
}


