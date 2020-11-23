#ifndef ELECTRICLIB_H
#define ELECTRICLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "TethysMathLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

class ElectroAnalysis{
private:
	std::ofstream data_electro;
	vector<float> TmpArr;
	vector<float> NetQ;
	vector<float> DipX;
	vector<float> DipY;
	vector<float> DipVarX;
	vector<float> DipVarY;
	vector<float> DipVarVarX;
	vector<float> DipVarVarY;
	vector<float> AvgCurDS;
	vector<float> VoltDS;
	vector<float> CurS;
	vector<float> CurD;
	vector<float> AvgCurHall;
	vector<float> PowOhm;
	vector<float> EngCap;
	vector<float> PowCap;

public:
	void CloseElectroFile();
	void CreateElectroFile(const GrapheneFluid1D& graphene);
	void CreateElectroFile(const string &infix_string);
	void WriteElectroFile(float t,const GrapheneFluid1D& graphene);
	void WriteElectroFile();
	static float NetCharge(const GrapheneFluid1D& graphene);
	static float NetCharge(const GrapheneFluid2D& graphene);
	static float OhmPower(const GrapheneFluid1D& graphene);
	static float OhmPower(const GrapheneFluid2D& graphene);
	static float AverageCurrent(const GrapheneFluid1D& graphene);
	static float AverageHallCurrent(const GrapheneFluid2D& graphene);
	static float AverageDirectCurrent(const GrapheneFluid2D& graphene);
	static float DrainCurrent(const GrapheneFluid2D& graphene);
	static float SourceCurrent(const GrapheneFluid2D& graphene);
	static float DrainToSourceVoltage(const GrapheneFluid2D& graphene);
	static float ElectricDipole(const GrapheneFluid1D& graphene);
	static float ElectricDipoleX(const GrapheneFluid2D& graphene);
	static float ElectricDipoleY(const GrapheneFluid2D& graphene);
	static float ElectricDipoleVariation(const GrapheneFluid1D& graphene);
	static float ElectricDipoleVariationX(const GrapheneFluid2D& graphene);
	static float ElectricDipoleVariationY(const GrapheneFluid2D& graphene);

	void ComputeElectroBase(float t, const GrapheneFluid2D &graphene);
	void ComputeElectroDerived();

	static void BannerDisplay(const GrapheneFluid2D &graphene);
};

#endif
