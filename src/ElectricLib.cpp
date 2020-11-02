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

void ElectroAnalysis::CreateElectroFile(const string &infix_string){
	std::string electrofile;
	electrofile = "electro_2D_" + infix_string + ".dat" ;
	data_electro.open (electrofile);
	data_electro << scientific;
}

void ElectroAnalysis::WriteElectroFile(float t,const GrapheneFluid1D& graphene){//TODO por em conformidade com o 2d talvez
	float q_net = this->NetCharge(graphene);
	float i_avg = this->AverageCurrent(graphene);
	float p_ohm = this->OhmPower(graphene);
	float dipole_var=this->ElectricDipoleVariation(graphene);
	float dipole=this->ElectricDipole(graphene);
	data_electro << t << "\t" << q_net << "\t" << i_avg << "\t" << q_net * q_net * 0.5 << "\t" << p_ohm << "\t" << dipole << "\t" << dipole_var << "\n";
}

void ElectroAnalysis::ComputeElectroBase(float t, const GrapheneFluid2D& graphene){
	TmpArr.push_back(t);
	NetQ.push_back(this->NetCharge(graphene));
	DipX.push_back(this->ElectricDipoleX(graphene));
	DipY.push_back(this->ElectricDipoleY(graphene));
	AvgCurDS.push_back(this->AverageDirectCurrent(graphene));
	AvgCurHall.push_back(this->AverageHallCurrent(graphene));
	PowOhm.push_back(this->OhmPower(graphene));
}


void ElectroAnalysis::ComputeElectroDerived() {
	float dt;
	dt=TmpArr.back()/DipX.size();
	EngCap.resize(DipX.size());
	PowCap.resize(DipX.size());
	transform(NetQ.begin(), NetQ.end(), EngCap.begin(), [](float &c){ return 0.5f*c*c; });
	Convolve_Gauss(1,5,1.0,EngCap.data(),PowCap.data(),EngCap.size());
	transform(PowCap.begin(), PowCap.end(), PowCap.begin(), [dt](float &c){ return c/dt; });
	DipVarX.resize(DipX.size());
	DipVarVarX.resize(DipX.size());
	DipVarY.resize(DipY.size());
	DipVarVarY.resize(DipY.size());
	Convolve_Gauss(1,5,1.0,DipX.data(),DipVarX.data(),DipX.size());
	Convolve_Gauss(1,5,1.0,DipVarX.data(),DipVarVarX.data(),DipX.size());
	Convolve_Gauss(1,5,1.0,DipY.data(),DipVarY.data(),DipY.size());
	Convolve_Gauss(1,5,1.0,DipVarY.data(),DipVarVarY.data(),DipY.size());
	transform(DipVarX.begin(), DipVarX.end(), DipVarX.begin(), [dt](float &c){ return c/dt; });
	transform(DipVarVarX.begin(), DipVarVarX.end(), DipVarVarX.begin(), [dt](float &c){ return c/(dt*dt); });
	transform(DipVarY.begin(), DipVarY.end(), DipVarY.begin(), [dt](float &c){ return c/dt; });
	transform(DipVarVarY.begin(), DipVarVarY.end(), DipVarVarY.begin(), [dt](float &c){ return c/(dt*dt); });
}


void ElectroAnalysis::WriteElectroFile() {
	try{
		if(NetQ.empty() || DipVarY.empty()){
			throw "Nothing to save on output file";
		}
		for(size_t i = 0; i < NetQ.size(); ++i) {
			data_electro << TmpArr[i] << "\t"
			             << NetQ[i] << "\t"
			             << AvgCurDS[i] << "\t"
			             << AvgCurHall[i] << "\t"
			             << PowOhm[i] << "\t"
			             << PowCap[i] << "\t"
			             << DipX[i] << "\t"
			             << DipVarX[i] << "\t"
			             << DipVarVarX[i] << "\t"
			             << DipY[i] << "\t"
			             << DipVarY[i] << "\t"
			             << DipVarVarY[i] << "\n";
		}
	}catch (const char* msg) {
		cerr << msg  <<"\nExiting"<< endl;
		exit(EXIT_FAILURE);
	}
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
		p += dx*static_cast<float>(2*j-2)*graphene.DenCor[2 * j - 2] + 4 * dx * static_cast<float>(2 * j - 1) * graphene.DenCor[2 * j - 1] + dx * static_cast<float>(2 * j) * graphene.DenCor[2 * j];
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
		auto i=static_cast<float>(divresult.rem);
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
		auto j=static_cast<float>(divresult.quot);
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
		pos=j*graphene.SizeX();
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}

float ElectroAnalysis::DrainCurrent(const GrapheneFluid2D& graphene) {
	int size1d = graphene.SizeY();
	float vector[size1d];
	int pos;
	for(int j=0;j<size1d;j++){
		pos=(j+1)*graphene.SizeX()-1;
		vector[j] = graphene.CurX[pos];
	}
	return Integral_1_D(graphene.SizeY(), graphene.GetDy(), vector);
}

void ElectroAnalysis::BannerDisplay(const GrapheneFluid2D &graphene) {
	graphene.BannerDisplay();
	cout<<"┌─────────────────────────────────────────────────────────────────────────┐\n";
	cout<<"│                \033[3mComputation of the Electronic Quantities\033[0m                 │\n";
	cout<<"└─────────────────────────────────────────────────────────────────────────┘\n";

	cout<<"Detected parameters: \n";
	cout << "Fermi velocity\t\033[1mvF\t" << graphene.GetVelFer() << " v\342\202\200\033[0m\n";
	cout << "Phase velocity\t\033[1mS'\t" << graphene.PhaseVel() << " v\342\202\200\033[0m\n";
	cout << "Viscosity \t\033[1m\316\267\t" << graphene.GetKinVis() << "\033[0m\n";
	if (graphene.GetKinVis() != 0.0) {
		cout << "Reynolds n. \t\033[1mRe\t" << 1.0 / graphene.GetKinVis() << "\033[0m\n";
	}
	cout << "Collision rate \t\033[1m\316\275\t" << graphene.GetColFreq()<< " v\342\202\200/L\033[0m\n";
	cout << "Cyclotron frequency \t\033[1m\317\211c\t" << graphene.GetCycFreq() << " v\342\202\200/L\n\033[0m\n";

}

