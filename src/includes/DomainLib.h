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
private:
	int size_x;
	int size_y;
	public:
		Domain(int Nx,int Ny);
		~Domain();
		bool *dom;///< Domain that contains the fluid
        void fill_Domain(function<float(float)> f_top,function<float(float)> f_bottom);///< fill the domain with the correct points
};
#endif