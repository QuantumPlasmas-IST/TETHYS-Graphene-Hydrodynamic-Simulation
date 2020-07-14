#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"

//TODO override the boundary conditions functions to Fluid1D and Fluid2D objects
class BoundaryCondition {
	public : 
	void XFree(GrapheneFluid1D& graphene);
	void XPeriodic(GrapheneFluid1D& graphene);
	void XFree(GrapheneFluid2D& graphene);
	void XPeriodic(GrapheneFluid2D& graphene);
	void YFree(GrapheneFluid2D& graphene);
	void YPeriodic(GrapheneFluid2D& graphene);
	void YClosedFreeSlip(GrapheneFluid2D& graphene);
	void YClosedNoSlip(GrapheneFluid2D& graphene);

	class DyakonovShur;
	class Dirichlet; 
};	

class  BoundaryCondition::DyakonovShur { 
	public:
	void X(GrapheneFluid1D& graphene);
	void X(GrapheneFluid2D& graphene);
	void YFree(GrapheneFluid2D& graphene);
	void YPeriodic(GrapheneFluid2D& graphene);
	void YClosedFreeSlip(GrapheneFluid2D& graphene);
	void YClosedNoSlip(GrapheneFluid2D& graphene);
};

class  BoundaryCondition::Dirichlet { 
	public: 	
	void Density(GrapheneFluid1D& graphene, float left, float right);
	void Density(GrapheneFluid2D& graphene, float left, float right, float top, float bottom);
	void VelocityX(GrapheneFluid1D& graphene, float left, float right);
	void MassFluxX(GrapheneFluid2D& graphene, float left, float right, float top, float bottom);
	void MassFluxY(GrapheneFluid2D& graphene, float left, float right, float top, float bottom);
};

#endif

