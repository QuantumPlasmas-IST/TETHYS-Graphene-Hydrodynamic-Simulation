#ifndef BOUNDARYLIB_H
#define BOUNDARYLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "Tethys2DLib.h"


class BoundaryCondition {
	public :
	void XFree(Fluid1D& fluid_class);
	void XFree(Fluid2D& fluid_class);
	void XPeriodic(Fluid1D& fluid_class);
	void XPeriodic(Fluid2D& fluid_class);
	void YFree(Fluid2D& fluid_class);
	void YPeriodic(Fluid2D& fluid_class);
	void YClosedFreeSlip(Fluid2D& fluid_class);
	void YClosedNoSlip(Fluid2D& fluid_class);

	class DyakonovShur;
	class Dirichlet; 
};	

class  BoundaryCondition::DyakonovShur { 
	public:
	void X(GrapheneFluid1D& fluid_class);
	void X(GrapheneFluid2D& fluid_class);
	void YFree(GrapheneFluid2D& fluid_class);
	void YPeriodic(GrapheneFluid2D& fluid_class);
	void YClosedFreeSlip(GrapheneFluid2D& fluid_class);
	void YClosedNoSlip(GrapheneFluid2D& fluid_class);
};

class  BoundaryCondition::Dirichlet { 
	public: 	
	void Density(Fluid1D& fluid_class, float left, float right);
	void Density(Fluid2D& fluid_class, float left, float right, float top, float bottom);
	void VelocityX(Fluid1D& fluid_class, float left, float right);
	void MassFluxX(Fluid2D& fluid_class, float left, float right, float top, float bottom);
	void MassFluxY(Fluid2D& fluid_class, float left, float right, float top, float bottom);
};

#endif

