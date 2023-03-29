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
	public:
        void condition_Edge(Fluid2D &fluid_class);///< The limit between   
};

#endif