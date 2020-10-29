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

	cout << "OPENING HDF5 file" <<endl;
	graph.OpenHdf5File(input_file_name);
	cout << "DONE" <<endl;

	elec.CreateElectroFile(graph);


	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		elec.WriteElectroFile(GrapheneFluid2D::TimeStamp,graph);
	}


	int numero_total=graph.GrpDen->getNumObjs();
	float time[numero_total];
	float dpx[numero_total];
	float ddpx[numero_total];
	float px[numero_total];
	for(hsize_t i=0; i < numero_total; i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		time[i]=GrapheneFluid2D::TimeStamp;

		px[i]=elec.ElectricDipoleX(graph);
		dpx[i]=0.0f;
		ddpx[i]=0.0f;
	}


	std::ofstream data_teste;
	std::string infix = graph.GetInfix();
	std::string testefile;
	testefile = "TESTE_2D_" + infix + ".dat" ;
	data_teste.open (testefile);
	data_teste << scientific;

Convolve_Gauss(1,5,1.0,px,dpx,numero_total);
Convolve_Gauss(1,5,1.0,dpx,ddpx,numero_total);
float diferenca_tempo;
	//diferenca_tempo=time[11]-time[10];
	diferenca_tempo=time[numero_total-1]/numero_total;
	for(hsize_t i=0; i < numero_total; i++){
		data_teste <<time[i]<<"\t"<< px[i] <<"\t"<< dpx[i]/diferenca_tempo <<"\t"<<dpx[i]/(diferenca_tempo*diferenca_tempo)<<endl;
	}

	return 0;
}

