/*!@file
 * @brief Header file for BC base class
 */

#ifndef EDGE_H
#define EDGE_H

#include <H5Cpp.h>
/*#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"
#include "BoundaryLib.h"*/
//#include "includes/DomainLib.h"
#include "DomainLib.h"

class Edge{
	private:
		int size_x;
		int size_y;
		
	public:
		Edge();
		Edge(int Nx, int Ny /*,bool * dom*/);
	    ~Edge();

		Domain D;

//		bool *edg;///< The limit between the domain and the outside
		std::vector <bool> edg;
		std::vector <int> edgint;
		
		void set_size_x(int);
		void set_size_y(int);
		void set_Edge();///< The limit between
};

#endif