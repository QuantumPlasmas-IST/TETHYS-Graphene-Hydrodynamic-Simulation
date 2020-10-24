#include "Tethys2DLib.h"
#include "ElectricLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

//const H5std_string FILE_NAME( "hdf5_2D_TEST.h5" );

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

	GrapheneFluid2D graph(parameters);

	cout << graph.SizeX()<<"\t"<<graph.SizeY()<<"\n";
return 0;
}