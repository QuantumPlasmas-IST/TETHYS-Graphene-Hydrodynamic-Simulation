#include "includes/BoundaryLib.h"
#include "includes/GeometryLib.h"

using namespace H5;
using namespace std;

Geometry::Geometry(int Nx, int Ny) : TethysBase{Nx,Ny,2} {
    size_x = Nx;
    size_y = Ny;
    fronteira.set_size_x(Nx);
    fronteira.set_size_y(Ny);
    dominio.set_size_x(Nx);
    dominio.set_size_y(Ny);
    fronteira.D.dom.resize(Nx*Ny);
    fronteira.edg.resize(Nx*Ny);
    dominio.dom.resize(Nx*Ny);
}

Geometry::~Geometry() = default;
/*
float Geometry::bottom_margin(float x){
    return (0.2*x);
}

float Geometry::top_margin(float x){
    return (0.3*x);
}
*/

void Geometry::SaveGeometry() {
//    vector <float> v;
//    v.resize(Nx*Ny);
    float w[200*400];
    cout << "cout" << endl;
    for(int k = 0; k < Nx*Ny; k++){
        cout << k << endl;
        if(fronteira.D.dom[k] == true){
            w[k] = 1;
        }else{
            if(fronteira.edg[k] == true){
                w[k] = 7;
            }else{
                w[k] = 0;
            }
        }
    }
    
/*	DataSet dataset_vel_snd = GrpDat->createDataSet("Sound velocity", HDF5FLOAT, *DataspaceVelSnd);
	dataset_vel_snd.write(w, HDF5FLOAT);
	dataset_vel_snd.close();*/
    cout << "oo" << endl;
    CreateHdf5File();          ///< creates the HDF5 files with the necessary structure
    cout << "debug 0" << endl;
  
    cout << "debug 1" << endl;

    DataSet dataset_vel_snd = GrpDat->createDataSet("Domain", HDF5FLOAT, *DataspaceVelSnd);
	cout << "debug 2" << endl;
    dataset_vel_snd.write(w, HDF5FLOAT);
	cout << "debug 3" << endl;
    dataset_vel_snd.close();
    cout << "debug 4" << endl;
    CloseHdf5File();  

}

/*bool Geometry::dom(int x, int y){
}

bool Geometry::edg(int x, int y){
}*/