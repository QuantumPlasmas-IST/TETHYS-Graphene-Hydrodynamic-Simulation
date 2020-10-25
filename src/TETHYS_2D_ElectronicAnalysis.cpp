#include "Tethys2DLib.h"
#include "ElectricLib.h"
#include "iomanip"
#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

int main(int argc, char **argv){
	std::cout << std::fixed;
	std::cout << std::setprecision(8);
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
	cout <<"graphene RANK "<<graph.Rank()<<endl;
	cout << "OPENING HDF5" <<endl;
	graph.OpenHdf5File("hdf5_2D_TEST.h5");
	cout << "DONE" <<endl;

	cout <<"Total number of datasets "<< graph.GrpDen->getNumObjs()<<endl;

	float t;
	//DataSet dataset_den;
	/*
	DataSpace* DataspaceDen;    // dataspace for EACH Density snapshots
	DataSpace* DataspaceVelX;   // dataspace for EACH Velocity X snapshots
	DataSpace* DataspaceVelY;   // dataspace for EACH Velocity Y snapshots
	*/

	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){

		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));


		cout << t<<"\t"<<graph.TimeStamp<<"\n";
		for(int k =0;k<5;k++){
			cout << graph.Den[k] <<"\t";
		}
		cout<<"\n";
		for(int k =0;k<5;k++){
			cout << graph.VelX[k] <<"\t";
		}
		cout<<"\n";
		for(int k =0;k<5;k++){
			cout << graph.VelY[k] <<"\t";
		}
		cout<<"\n\n";


	}


return 0;
}

