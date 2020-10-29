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
public:
	vector<float> TmpArr;
	vector<float> NetQ;
	vector<float> DipX;
	vector<float> DipY;
	vector<float> DipVarX;
	vector<float> DipVarY;
	vector<float> DipVarVarX;
	vector<float> DipVarVarY;
	vector<float> AvgCurDS;
	vector<float> AvgCurHall;
	vector<float> PowOhm;

	void CloseElectroFile();
	void CreateElectroFile(const GrapheneFluid1D& graphene);
	void CreateElectroFile(const GrapheneFluid2D& graphene);
	void WriteElectroFile(float t,const GrapheneFluid1D& graphene);
	void WriteElectroFile(float t,const GrapheneFluid2D& graphene);
	float NetCharge(const GrapheneFluid1D& graphene);
	float NetCharge(const GrapheneFluid2D& graphene);
	float OhmPower(const GrapheneFluid1D& graphene);
	float OhmPower(const GrapheneFluid2D& graphene);
	float AverageCurrent(const GrapheneFluid1D& graphene);
	float AverageHallCurrent(const GrapheneFluid2D& graphene);
	float AverageDirectCurrent(const GrapheneFluid2D& graphene);
	float DrainCurrent(const GrapheneFluid2D& graphene);
	float SourceCurrent(const GrapheneFluid2D& graphene);
	float DrainToSourceVoltage(const GrapheneFluid2D& graphene);
	float ElectricDipole(const GrapheneFluid1D& graphene);
	float ElectricDipoleX(const GrapheneFluid2D& graphene);
	float ElectricDipoleY(const GrapheneFluid2D& graphene);
	float ElectricDipoleVariation(const GrapheneFluid1D& graphene);
	float ElectricDipoleVariationX(const GrapheneFluid2D& graphene);
	float ElectricDipoleVariationY(const GrapheneFluid2D& graphene);

	void ComputeElectro(float t, const GrapheneFluid2D &graphene);
};

#endif
