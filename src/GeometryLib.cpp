#include "includes/BoundaryLib.h"
#include "includes/GeometryLib.h"

float Geometry::bottom_margin(float x){
    return (0.2*x);
}

float Geometry::top_margin(float x){
    return (0.3*x);
}

void Geometry::SaveGeometry() {

	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(dominio.dom, HDF5FLOAT);
	dataset_vel_snd.close();
}

/*bool Geometry::dom(int x, int y){
}

bool Geometry::edg(int x, int y){
}*/