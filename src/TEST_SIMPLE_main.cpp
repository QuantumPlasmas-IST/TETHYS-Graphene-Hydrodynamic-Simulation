#include "includes/Fluid1DLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid1DLib.h"
#include "includes/InitialConditionLib.h"
#include "includes/FeedbackBoundaryLib.h"
#include "BoundaryLib.h"
#include "GeometryLib.h"

#include <functional>

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

float ff_top(float x){
    return 0.2*x;
}

float ff_bottom(float x){
    return 0.3*x;
}

int main(){
    int NX = 400;
    int NY = 200;
    cout << "w" << endl;
    ofstream myfile;
    myfile.open("wow.txt");
    Geometry Geo(NX,NY);
    cout << "w" << endl;
    Geo.fronteira.D.set_Domain(ff_top,ff_bottom);
    cout << "w" << endl;
    Geo.fronteira.set_Edge();
    cout << "w" << endl;
    for(int k = 0; k < NX*NY-1; k++){
//        cout << Geo.fronteira.D.dom[k];
        if(Geo.fronteira.edg[k] == false){
            myfile << Geo.fronteira.D.dom[k];
        }
        if(Geo.fronteira.edg[k] == true){
            myfile << "L"; 
        }
        if((k+1)%NX == 0){
//            cout << " k = " << k << endl;
            myfile << endl;
        }
    }
    cout << "writing" << endl;
    myfile.close();
    cout << "writing" << endl;
    Geo.SaveGeometry();
    cout << "writing" << endl;
//    textfile *meh;
//    meh.write();
    return 0;
}