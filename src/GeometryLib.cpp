#include "includes/BoundaryLib.h"
#include "includes/GeometryLib.h"

float Geometry::bottom_margin(float x){
    return (0.2*x);
}

float Geometry::top_margin(float x){
    return (0.3*x);
}

/*bool Geometry::dom(int x, int y){
}

bool Geometry::edg(int x, int y){
}*/