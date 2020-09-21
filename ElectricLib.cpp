#include "ElectricLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

using namespace H5;
using namespace std;


void ElectroAnalysis::CreateElectroFile(GrapheneFluid1D& graphene){
	std::string infix = graphene.GetInfix();
	std::string electrofile;
	electrofile = "electro_1D_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::CreateElectroFile(GrapheneFluid2D& graphene){
	std::string infix = graphene.GetInfix();
	std::string electrofile;
	electrofile = "electro_2D_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::WriteElectroFile(float t, GrapheneFluid1D& graphene){
	float q_net = this->NetCharge(graphene);
	float i_avg = this->AverageCurrent(graphene);
	float p_ohm = this->OhmPower(graphene);
	float dipole_var=this->ElectricDipoleVariation(graphene);
	float dipole=this->ElectricDipole(graphene);
	data_electro << t << "\t" << q_net << "\t" << i_avg << "\t" << q_net * q_net * 0.5 << "\t" << p_ohm << "\t" << dipole << "\t" << dipole_var << "\n";
}

void ElectroAnalysis::WriteElectroFile(float t, GrapheneFluid2D& graphene){
	data_electro << t << "\t" << "\n";
}

float ElectroAnalysis::NetCharge(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.DenCor);
}
float ElectroAnalysis::NetCharge(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.Den);
}

float ElectroAnalysis::AverageCurrent(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.CurCor);
}

float ElectroAnalysis::AverageCurrent(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.CurX);
}

float ElectroAnalysis::ElectricDipoleVariation(const GrapheneFluid1D& graphene){
	return Integral_1_D(graphene.SizeX(), graphene.GetDx(), graphene.CurCor);
}

float ElectroAnalysis::ElectricDipole(const GrapheneFluid1D& graphene){
	float p=0.0;
	float dx=graphene.GetDx();
	for(int j=1;j<graphene.SizeX()/2;j++){
		p += dx*(2*j-2)*graphene.DenCor[2 * j - 2] + 4 * dx * (2 * j - 1) * graphene.DenCor[2 * j - 1] + dx * (2 * j) * graphene.DenCor[2 * j];
	}
	p = p*graphene.GetDx()/3.0f;
	return p;
}

float ElectroAnalysis::OhmPower(const GrapheneFluid1D& graphene){
	float itg=0.0;
	for(int j=1;j<graphene.SizeX()/2;j++){
		itg += graphene.CurCor[2 * j - 2] * graphene.VelCor[2 * j - 2] + 4 * graphene.CurCor[2 * j - 1] * graphene.VelCor[2 * j - 1] + graphene.CurCor[2 * j] * graphene.VelCor[2 * j];
	}
	itg = itg*graphene.GetDx()/3.0f;
	return itg;
}

