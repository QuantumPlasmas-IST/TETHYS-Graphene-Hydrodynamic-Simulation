#include "Tethys2DLib.h"
#include "ElectricLib.h"



using namespace std;

int main(int argc, char **argv){

	string input_file_name;
	if(argc==2){
		input_file_name = argv[1];
	}else{
		cout <<"Input HDF5 file to read:\t";
		cin >> input_file_name;
	}
	cout<<"Reading file:\t"<< input_file_name <<endl;

	SetUpParameters parameters;
	parameters.ParametersFromHdf5File(input_file_name);

	cout <<"Sound velocity\t"<< parameters.SoundVelocity << endl;
	cout <<"Fermi velocity\t"<< parameters.FermiVelocity << endl;
	cout <<"Collision frequency\t"<< parameters.CollisionFrequency << endl;
	cout <<"Viscosity\t"<< parameters.ShearViscosity<< endl;
	cout <<"Cyclotron frequency\t"<< parameters.CyclotronFrequency << endl;
	cout <<"Aspect Ratio\t"<< parameters.AspectRatio << endl;
	cout <<"Dimensions\t"<< parameters.Length <<" x "<< parameters.Width << endl;
	cout <<"Grid\t"<< parameters.SizeX <<" x "<< parameters.SizeY << endl;

	GrapheneFluid2D graph(parameters);
	ElectroAnalysis elec;

	cout << "Importing HDF5 file...";
	graph.OpenHdf5File(input_file_name);
	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		elec.ComputeElectroBase(GrapheneFluid2D::TimeStamp, graph);
	}
	cout << " DONE" <<endl;
	elec.ComputeElectroDerived();
	elec.CreateElectroFile(graph);
	cout << "Writing outpu file...";
	elec.WriteElectroFile(GrapheneFluid2D::TimeStamp,graph);
	cout << " DONE" <<endl;

	return 0;
}

