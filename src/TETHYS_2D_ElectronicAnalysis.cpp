/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/Fluid2DLib.h"
#include "includes/ElectricLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/GrapheneFluid2DLib.h"


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

	ElectroAnalysis::BannerDisplay(graph);

	graph.OpenHdf5File(input_file_name);
	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		elec.ComputeElectroBase(GrapheneFluid2D::TimeStamp, graph);
	}
	elec.ComputeElectroDerived();
	elec.CreateElectroFile(graph.GetInfix());
	elec.WriteElectroFile();

	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}

