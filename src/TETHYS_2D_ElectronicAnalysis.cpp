#include "Tethys2DLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

//const H5std_string FILE_NAME( "hdf5_2D_TEST.h5" );
herr_t
file_info(hid_t loc_id, const char *name, void *opdata)
{
	/*
	 * Display group name.
	 */
	cout << "Name : " << name<< endl;

	return 0;
}

int main(int argc, char **argv){

	SetUpParameters parameters;
	parameters.ParametersFromHdf5File("hdf5_2D_TEST.h5");
	cout <<"NX\t"<< parameters.SizeX << endl;
	cout <<"NY\t"<< parameters.SizeY << endl;
	cout <<"Sound velocity\t"<< parameters.SoundVelocity << endl;
	cout <<"Fermi velocity\t"<< parameters.FermiVelocity << endl;
	cout <<"Collision frequency\t"<< parameters.CollisionFrequency << endl;
	cout <<"Viscosity\t"<< parameters.ShearViscosity<< endl;
	cout <<"Cyclotron frequency\t"<< parameters.CyclotronFrequency << endl;
	cout <<"Aspect Ratio\t"<< parameters.AspectRatio << endl;

	//GrapheneFluid2D graph(parameters);

	//cout << graph.SizeX()<<"\t"<<graph.SizeY()<<"\n";

	cout << "OPENING HDF5" <<endl;

	H5std_string  FILE_NAME( "hdf5_2D_TEST.h5" );
	H5File* hdf5_file;
	Group* grp_dat;
	Group* grp_den;
	hdf5_file = new H5File(FILE_NAME, H5F_ACC_RDONLY );
	grp_dat = new Group(hdf5_file->openGroup("/Data" ));
	grp_den = new Group(hdf5_file->openGroup("/Data/Density" ));

	cout <<"numero de coisos"<< grp_den->getNumObjs()<<endl;

	for(int i=0;i<grp_den->getNumObjs();i++){

	cout <<"Dataset ID:"<< i <<"\tName:"<<grp_den->getObjnameByIdx(i)<<endl;
	}


return 0;
}

