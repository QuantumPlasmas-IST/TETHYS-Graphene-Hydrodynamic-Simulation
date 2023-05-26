/*!@file
 * @brief Header file for BC base class
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <functional>
#include <H5Cpp.h>
#include "includes/TethysBaseLib.h"
#include "includes/TethysMathLib.h"
/*#include "includes/Fluid1DLib.h"
#include "includes/Fluid2DLib.h"
#include "includes/DiracGraphene2DLib.h"*/
#include "includes/DomainLib.h"
#include "includes/EdgeLib.h"
//include "BoundaryLib.h"

using namespace std;

class Geometry : public TethysBase{
        private:
                int size_x;
                int size_y;
        public:
                Geometry(int, int);
                ~Geometry();
                
//                float bottom_margin(float x);///< Function that describes the edge on the botton side
//                float top_margin(float x);///< Function that describes the edge on the botton side
//                bool dom(int x, int y);///< Domain that contains the fluid
//                bool edg(int x, int y);///< The limit between the domain and the outside  
		Domain dominio;
		Edge fronteira;

                void SaveGeometry();

};

#endif