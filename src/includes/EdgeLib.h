/*!@file
 * @brief Header file for BC base class
 */

#ifndef EDGE_H
#define EDGE_H

#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"
#include "BoundaryLib.h"
#include "includes/GeometryLib.h"

class Edge : public Geometry{
	private:
		int size_x;
		int size_y;
	public:
		Edge(int Nx, int Ny,bool * dom);
	    ~Edge();

		bool *edg;///< The limit between the domain and the outside

		void condition_Edge(function<float(float)> f_top,function<float(float)> f_bottom);///< The limit between
};

#endif