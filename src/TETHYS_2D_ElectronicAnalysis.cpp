#include "Tethys2DLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

const H5std_string FILE_NAME( "hdf5_2D_TEST.h5" );

int main(int argc, char **argv){

	H5File* file;
	Group* group;
	file = new H5File( FILE_NAME, H5F_ACC_RDONLY );
	group = new Group(file->openGroup("Data"));
	Attribute *attr_n_x = new Attribute(group->openAttribute("Number of spatial points x"));
	Attribute *attr_n_y = new Attribute(group->openAttribute("Number of spatial points y"));
	Attribute *attr_snd = new Attribute(group->openAttribute("Sound velocity"));
	Attribute *attr_fer = new Attribute(group->openAttribute("Fermi velocity"));
	Attribute *attr_vis = new Attribute(group->openAttribute("Kinetic viscosity"));
	Attribute *attr_cyc = new Attribute(group->openAttribute("Cyclotron frequency"));
	Attribute *attr_col = new Attribute(group->openAttribute("Collision frequency"));

	float test_float;
	int test_int;


	int npoints_x,npoints_y;
	float sound, fermi, coll, visco, cyclo;
	attr_n_x->read(attr_n_x->getDataType(), &npoints_x);
	attr_n_y->read(attr_n_y->getDataType(), &npoints_y);
	attr_snd->read(attr_snd->getDataType(), &sound);
	attr_fer->read(attr_fer->getDataType(), &fermi);
	attr_vis->read(attr_vis->getDataType(), &visco);
	attr_cyc->read(attr_cyc->getDataType(), &cyclo);
	attr_col->read(attr_col->getDataType(), &coll);

	cout <<"Sound velocity\t"<< sound << endl;
	cout <<"Fermi velocity\t"<< fermi << endl;
	cout <<"Collision frequency\t"<< coll << endl;
	cout <<"Viscosity\t"<< visco << endl;
	cout <<"Cyclotron frequency\t"<< cyclo << endl;
	cout <<"Aspect Ratio\t"<< (npoints_x-1)/(npoints_y-1) << endl;

	SetUpInput parameters(sound, fermi, coll, visco, cyclo, 1, (npoints_x-1)/(npoints_y-1) );
	GrapheneFluid2D	graph(npoints_x, npoints_y, parameters);

	delete group;
	delete file;
return 0;
}