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

	SetUpParameters parameters;
	parameters.ParametersFromHdf5File(input_file_name);

	GrapheneFluid2D graph(parameters);
	ElectroAnalysis elec;

	elec.BannerDisplay(graph);


	cout << "Importing HDF5 file...";
	graph.OpenHdf5File(input_file_name);
	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		elec.ComputeElectroBase(GrapheneFluid2D::TimeStamp, graph);
	}
	cout << " DONE" <<endl;
	elec.ComputeElectroDerived();
	elec.CreateElectroFile(graph.GetInfix());
	cout << "Writing outpu file...";
	elec.WriteElectroFile();
	cout << " DONE" <<endl;

	return 0;
}

