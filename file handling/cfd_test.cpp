
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
const H5std_string      DATASET_NAME( "Density" );
const int       NX = 100;                    // dataset dimensions
const int       NY = 100;
const int       RANK = 2;

int main (void)
{
	/*
	* Data initialization.
	*/
	int i, j;
	float data[NX][NY];          // buffer for data to write
	for (j = 0; j < NX; j++)
	{
	  for (i = 0; i < NY; i++)
		 data[j][i] = sin(i*0.1) + 2.0*sin(j*0.1);
	}
  
	/*
	 * Create a file.
	 */
	H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	/*
	* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	hsize_t     dimsf[2];              // dataset dimensions
	dimsf[0] = NX;
	dimsf[1] = NY;
	DataSpace dataspace( RANK, dimsf );

FloatType datatype(PredType::NATIVE_FLOAT);
DataSet dataset = file->createDataSet( DATASET_NAME, datatype, dataspace );
dataset.write( data, PredType::NATIVE_FLOAT );

	return 0;  // successfully terminated
}

