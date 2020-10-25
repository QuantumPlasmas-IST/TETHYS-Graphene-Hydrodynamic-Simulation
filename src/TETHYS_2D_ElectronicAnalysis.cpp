#include "Tethys2DLib.h"
#include "ElectricLib.h"



using namespace std;

int main(int argc, char **argv){
	//std::cout << std::fixed;
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
	ElectroAnalysis elec;
	cout <<"graphene RANK "<<graph.Rank()<<endl;
	cout << "OPENING HDF5" <<endl;
	graph.OpenHdf5File("hdf5_2D_TEST.h5");
	cout << "DONE" <<endl;

	cout <<"Total number of datasets "<< graph.GrpDen->getNumObjs()<<endl;



	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){

		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		cout <<graph.TimeStamp<<"\t"
		<< elec.NetCharge(graph)<<"\t"
		<< elec.AverageDirectCurrent(graph) <<"\t"
		<< elec.AverageHallCurrent(graph) <<"\t"
		<< elec.OhmPower(graph) <<"\t"
		<< elec.ElectricDipoleVariation(graph) <<"\t"
		<< elec.ElectricDipole(graph) <<"\n";
	}





	return 0;
}

