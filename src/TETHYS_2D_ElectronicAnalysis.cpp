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


	string infix = graph.GetInfix();
	string electrofile;
	electrofile = "RADIATION_2D_" + infix + ".dat" ;
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific;

	for(hsize_t i=0; i < graph.GrpDen->getNumObjs(); i++){
		graph.ReadSnapShot(graph.GrpDen->getObjnameByIdx(i));
		graph.VelocityToCurrent();
		data_electro << graph.TimeStamp << "\t"
		             << elec.NetCharge(graph) << "\t"
		             << elec.AverageDirectCurrent(graph) << "\t"
		             << elec.AverageHallCurrent(graph) << "\t"
		             << elec.OhmPower(graph) << "\t"
		             << elec.ElectricDipoleVariationX(graph) << "\t"
		             << elec.ElectricDipoleVariationY(graph) << "\t"
		             << elec.ElectricDipoleX(graph) << "\t"
					 << elec.ElectricDipoleY(graph) << "\n";
	}


	cout << "Benchmark Electrical functions"<<endl;
	for(int k=0;k<graph.SizeX()*graph.SizeY();k++){
		div_t divresult;
		divresult = div (k,graph.SizeX());
		int j=divresult.quot;
		int i=divresult.rem;
		float x = i*graph.GetDx();
		float y = j*graph.GetDy();
		graph.Den[k]  =  1.0f+sin(3.0f*MAT_PI*x)-5.0f*y*(y-1.0f)+7.0f*x;
		graph.VelX[k] =  pow(y*(y-1.0f),2);
		graph.VelY[k] =  0.2f*x;
	}
	graph.VelocityToCurrent();



	float net_charge = 16.0f/3.0f+2.0f/(3.0f*MAT_PI);
	float i_ds = 13.0f/70.0f+1.0f/(45.0f*MAT_PI);
	float i_hall = 39.0f/60.0f+4.0f/(60.0f*MAT_PI);
	float p_ohm = 233333.0f/307125.0f-296.0f/(2025.0f*MAT_PI*MAT_PI*MAT_PI)-1.0f/(900.0f*MAT_PI*MAT_PI)+12863.0f/(51975.0f*MAT_PI);
	float px = 7.0f/12.0f;
	float py = 0.0f;




	cout <<"Sim. dimensions  "<<graph.GetLengthX()<<" x "<< graph.GetLengthY()<<"\n";
	cout <<"Sim. points  "<<graph.SizeX()<<" x "<< graph.SizeY()<<"\n";
	std::cout << std::scientific;
	std::cout << std::setprecision(8);
	cout << "Quantity\t\t"<<"Numerical"<<"\t"<< "Analytical"<<"\t"<<"%error"<<"\n";
	cout << "Total charge\t\t"<<elec.NetCharge(graph)<<"\t"<< net_charge <<"\t"<< 100.0f*abs(elec.NetCharge(graph)-net_charge)/net_charge<<"\n";
	cout << "Avg. current\t\t"<<elec.AverageDirectCurrent(graph)<<"\t"<< i_ds <<"\t"<<  100.0f*abs(elec.AverageDirectCurrent(graph)-i_ds)/i_ds <<"\n";
	cout << "Avg. Hall current\t"<<elec.AverageHallCurrent(graph)<<"\t"<<  i_hall <<"\t"<< 100.0f*abs(elec.AverageHallCurrent(graph)-i_hall)/i_hall<<"\n";
	cout << "Ohm Power\t\t"<<elec.OhmPower(graph)<<"\t"<< p_ohm <<"\t"<< 100.0f*abs(elec.OhmPower(graph)-p_ohm)/p_ohm <<"\n";
	cout << "Elec. Dipole x\t\t"<<elec.ElectricDipoleX(graph)<<"\t"<< px <<"\t"<< 100.0f*abs(elec.ElectricDipoleX(graph)-px)/px <<"\n";
	cout << "Elec. Dipole y\t\t"<<elec.ElectricDipoleY(graph)<<"\t"<< py <<"\t"<< 100.0f*abs(elec.ElectricDipoleY(graph)-py)/py <<"\n";
	cout << "Elec. Dipole var x\t"<<elec.ElectricDipoleVariationX(graph)<<"\t"<<i_ds <<"\t"<<  100.0f*abs(elec.AverageDirectCurrent(graph)-i_ds)/i_ds <<"\n";
	cout << "Elec. Dipole var y\t"<<elec.ElectricDipoleVariationY(graph)<<"\t"<<i_hall <<"\t"<< 100.0f*abs(elec.AverageHallCurrent(graph)-i_hall)/i_hall <<"\n";


	return 0;
}

