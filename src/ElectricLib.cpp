#include "ElectricLib.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

using namespace H5;
using namespace std;


void ElectroAnalysis::CreateElectroFile(const GrapheneFluid1D& graphene){
	std::string infix = graphene.GetInfix();
	std::string electrofile;
	electrofile = "electro_1D_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::CreateElectroFile(const GrapheneFluid2D& graphene){
	std::string infix = graphene.GetInfix();
	std::string electrofile;
	electrofile = "electro_2D_" + infix + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::WriteElectroFile(float t,const GrapheneFluid1D& graphene){
	float q_net = this->NetCharge(graphene);
	float i_avg = this->AverageCurrent(graphene);
	float p_ohm = this->OhmPower(graphene);
	float dipole_var=this->ElectricDipoleVariation(graphene);
	float dipole=this->ElectricDipole(graphene);
	data_electro << t << "\t" << q_net << "\t" << i_avg << "\t" << q_net * q_net * 0.5 << "\t" << p_ohm << "\t" << dipole << "\t" << dipole_var << "\n";
}

void ElectroAnalysis::WriteElectroFile(float t,const GrapheneFluid2D& graphene){
	float q_net = this->NetCharge(graphene);
	float i_ds   = this->AverageDirectCurrent(graphene);
	float i_hall = this->AverageHallCurrent(graphene);
	float p_ohm = this->OhmPower(graphene);
	float dipole_var_x= this->ElectricDipoleVariationX(graphene);
	float dipole_x= this->ElectricDipoleX(graphene);
	float dipole_var_y= this->ElectricDipoleVariationY(graphene);
	float dipole_y= this->ElectricDipoleY(graphene);
	data_electro << t << "\t"<<q_net << "\t"<<i_ds<<"\t"<<i_hall<< "\t" << p_ohm << "\t" << dipole_x << "\t" << dipole_var_x << dipole_y << "\t" << dipole_var_y << "\n";
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

float ElectroAnalysis::AverageDirectCurrent(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.CurX);
}

float ElectroAnalysis::AverageHallCurrent(const GrapheneFluid2D& graphene){
	return Integral_2_D(graphene.SizeX(),graphene.SizeY(), graphene.GetDx(),graphene.GetDy(), graphene.CurY);
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

float ElectroAnalysis::OhmPower(const GrapheneFluid2D& graphene){
	int size = graphene.SizeX()*graphene.SizeX();
	float square_current_density[size];
	float jx,jy;
	for(int c=0;c<size;c++){
		jx=graphene.CurX[c];
		jy=graphene.CurY[c];
		square_current_density[c] = jx*jx+jy*jy;
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), square_current_density);
}

float ElectroAnalysis::ElectricDipoleX(const GrapheneFluid2D &graphene) {
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	float rx;
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		int i=divresult.rem;
		rx = i*graphene.GetDx()-0.5f*graphene.GetLengthX();
		vector[c] = rx*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}
float ElectroAnalysis::ElectricDipoleY(const GrapheneFluid2D &graphene) {
	int size = graphene.SizeX()*graphene.SizeY();
	float vector[size];
	float ry;
	for(int c=0;c<size;c++){
		div_t divresult;
		divresult = div (c,graphene.SizeX());
		int j=divresult.quot;
		ry = j*graphene.GetDy()-0.5f*graphene.GetLengthY();
		vector[c] = ry*graphene.Den[c];
	}
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), vector);
}
float ElectroAnalysis::ElectricDipoleVariationX(const GrapheneFluid2D &graphene) {
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), graphene.CurX );
}
float ElectroAnalysis::ElectricDipoleVariationY(const GrapheneFluid2D &graphene) {
	return Integral_2_D(graphene.SizeX(), graphene.SizeY(), graphene.GetDx(), graphene.GetDy(), graphene.CurY );
}

void ElectroAnalysis::CloseElectroFile() {
	data_electro.close();
}

float ElectroAnalysis::DrainToSourceVoltage(const GrapheneFluid2D& graphene) {
	return 0;
}

float ElectroAnalysis::SourceCurrent(const GrapheneFluid2D& graphene) {
	int size1d = graphene.SizeY();
	float vector[size1d];
	int pos;
	for(int j=0;j<size1d;j++){
		pos=graphene.SizeX()-1+j*graphene.SizeX();
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}

float ElectroAnalysis::DrainCurrent(const GrapheneFluid2D& graphene) {
	int size1d = graphene.SizeY();
	float vector[size1d];
	int pos;
	for(int j=0;j<size1d;j++){
		pos=j*graphene.SizeX();
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}
// i+j*graphene.SizeY()