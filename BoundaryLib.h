#ifndef BONDLIB_H
#define BONDLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"


class BoundaryCondition {
	public : 
	void XFree(GrapheneFluid1D& graphene);
	void XPeriodic(GrapheneFluid1D& graphene);
	void XFree(GrapheneFluid2D& graphene);
	void XPeriodic(GrapheneFluid2D& graphene);
	void YFree(GrapheneFluid2D& graphene);
	void YPeriodic(GrapheneFluid2D& graphene);

	class DyakonovShur;
	class Dirichlet; 
};	

class  BoundaryCondition::DyakonovShur { 
	public:
	void X(GrapheneFluid1D& graphene);
	void X(GrapheneFluid2D& graphene);
	void YFree(GrapheneFluid2D& graphene);
	void YPeriodic(GrapheneFluid2D& graphene);
};

class  BoundaryCondition::Dirichlet { 
	public: 	
	void Density(GrapheneFluid1D& graphene, float L, float R);
	void Density(GrapheneFluid2D& graphene, float L, float R, float T, float B);
	void VelocityX(GrapheneFluid1D& graphene, float L, float R);
	void VelocityX(GrapheneFluid2D& graphene, float L, float R, float T, float B);
	void VelocityY(GrapheneFluid2D& graphene, float L, float R, float T, float B);
};

#endif
