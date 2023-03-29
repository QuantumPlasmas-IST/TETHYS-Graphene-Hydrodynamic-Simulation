/*!@file
 * @brief Header file for BC base class
 */

#ifndef DOMAIN_H
#define DOMAIN_H


#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"
#include "includes/BoundaryLib.h"
#include "includes/GeometryLib.h"

class Domain : public Geometry{
	public:
        void fill_Domain(Fluid2D &fluid_class);///< fill the domain with the correct points   
};
#endif