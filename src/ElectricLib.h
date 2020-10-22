#ifndef ELECTRICLIB_H
#define ELECTRICLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

class ElectroAnalysis{
private:
	std::ofstream data_electro;
public:
	void CreateElectroFile(const GrapheneFluid1D& graphene);
	void CreateElectroFile(const GrapheneFluid2D& graphene);
	void WriteElectroFile(float t,const GrapheneFluid1D& graphene);
	void WriteElectroFile(float t,const GrapheneFluid2D& graphene);
	float NetCharge(const GrapheneFluid1D& graphene);
	float NetCharge(const GrapheneFluid2D& graphene);
	float OhmPower(const GrapheneFluid1D& graphene);
	//float OhmPower(GrapheneFluid2D& graphene);
	float AverageCurrent(const GrapheneFluid1D& graphene);
	float AverageHallCurrent(const GrapheneFluid2D& graphene);
	float AverageDirectCurrent(const GrapheneFluid2D& graphene);
	float ElectricDipole(const GrapheneFluid1D& graphene);
	//float ElectricDipole(GrapheneFluid2D& graphene);
	float ElectricDipoleVariation(const GrapheneFluid1D& graphene);
	//float ElectricDipoleVariation(GrapheneFluid2D& graphene);
};

#endif
