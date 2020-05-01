
/*
 *  This example writes a dataset to a new HDF5 file.
*/


#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>

#include "H5Cpp.h"

#include "dyakonovshur.h"

using namespace H5;
using namespace std;

const H5std_string      FILE_NAME( "CFD.h5" );
const H5std_string      DATASET_DEN( "Density" );
const H5std_string      DATASET_VEL( "Velocity" );
const int       NX = 200;                    // dataset dimensions
const int       NY = 200;
const int       RANK = 1;

int main (void)
{
	/*
	* Data initialization.
	*/
	
	float data_den[NX];          // buffer for data to write
	float data_vel[NX];          // buffer for data to write
	for (int j = 0; j < NX; j++)
	{
	  for (int i = 0; i < NY; i++){
		 data_den[j] = sin(i*0.1) + 2.0*sin(j*0.1);
		 data_vel[j] = i*0.1-j*0.1;
		}
	}
  
	/*
	 * Create a file.
	 */
	H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	/*
	 * Create a group in the file
	 */
	Group* group_data = new Group( file->createGroup( "/Data" ));
	Group* group_den = new Group( file->createGroup( "/Data/Density" ));
	Group* group_vel = new Group( file->createGroup( "/Data/Velocity" ));




	/*
	 * Create attributes 
	 */	
	float  vel_snd, vel_fer , col_freq;
	cout << "Define S value: ";
	cin >> vel_snd;
	cout << "Define vF value: ";
	cin >> vel_fer;
	cout << "Define collision frequency: ";
	cin >> col_freq; 
	 
	//const int	DIM1 = 1;
    hsize_t dims[1] = { 1 };
	// Create the data space for the attribute.
	DataSpace attr_dataspace = DataSpace (1, dims );
	// Create a group attribute. 
	Attribute attribute1 = group_data->createAttribute( "S parameter", PredType::NATIVE_FLOAT, attr_dataspace);
	Attribute attribute2 = group_data->createAttribute( "Fermi velocity", PredType::NATIVE_FLOAT, attr_dataspace);
	Attribute attribute3 = group_data->createAttribute( "Collision frequency", PredType::NATIVE_FLOAT, attr_dataspace);
	// Write the attribute data. 
	attribute1.write( PredType::NATIVE_FLOAT, &vel_snd);
	attribute2.write( PredType::NATIVE_FLOAT, &vel_fer);
	attribute3.write( PredType::NATIVE_FLOAT, &col_freq);
	
	
	/*
	* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	hsize_t     dimsf[1];              // dataset dimensions
	dimsf[0] = NX;
	//dimsf[1] = NY;
	DataSpace dataspace_den( RANK, dimsf );
	DataSpace dataspace_vel( RANK, dimsf );


	FloatType datatype_flt(PredType::NATIVE_FLOAT);
	
	DataSet dataset_den = group_den->createDataSet(DATASET_DEN, PredType::NATIVE_FLOAT, dataspace_den );
	dataset_den.write( data_den, PredType::NATIVE_FLOAT );



//	DataSet dataset_vel = group_vel->createDataSet(DATASET_VEL , datatype_flt, dataspace_vel );
//	dataset_vel.write( data_vel, PredType::NATIVE_FLOAT );


	for(int t=0 ;t<=5; t++){
		
		for (int j = 0; j < NX; j++)
		{
			for (int i = 0; i < NY; i++){
				data_vel[j] = i*0.1-j*0.1*t;
			}
		}
		
		string str_time = to_string(t);
		//string name_dateset = "/Data/"+str_time;
		
		DataSet dataset_vel = group_vel->createDataSet( str_time , datatype_flt, dataspace_vel );
		dataset_vel.write( data_vel, PredType::NATIVE_FLOAT );

		
	
	
	}
	
	
    delete &dataset_vel
    delete file;
	return 0;  // successfully terminated
}

