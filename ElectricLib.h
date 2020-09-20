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
	void CreateElectroFile(GrapheneFluid1D& graphene);
	void CreateElectroFile(GrapheneFluid2D& graphene);
	void WriteElectroFile(float t,GrapheneFluid1D& graphene);
	void WriteElectroFile(float t,GrapheneFluid2D& graphene);
	float NetCharge(GrapheneFluid1D& graphene);
	float NetCharge(GrapheneFluid2D& graphene);
	float OhmPower(GrapheneFluid1D& graphene);
	//float OhmPower(GrapheneFluid2D& graphene);
	float AverageCurrent(GrapheneFluid1D& graphene);
	//float AverageCurrent(GrapheneFluid2D& graphene);
	float ElectricDipole(GrapheneFluid1D& graphene);
	//float ElectricDipole(GrapheneFluid2D& graphene);
	float ElectricDipoleVariation(GrapheneFluid1D& graphene);
	//float ElectricDipoleVariation(GrapheneFluid2D& graphene);
};

#endif
