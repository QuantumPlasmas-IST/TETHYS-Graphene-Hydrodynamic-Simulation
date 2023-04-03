#include "includes/DomainLib.h"

Domain::Domain(){}

Domain::Domain(int Nx,int Ny /*, function<float(float)> f_top,function<float(float)> f_bottom*/){
    size_x = Nx;
    size_y = Ny;
    dom.resize(Nx*Ny);
}

Domain::~Domain() = default;

void Domain::set_size_x(int Nx){
    size_x = Nx;
}

void Domain::set_size_y(int Ny){
    size_y = Ny;
}

void Domain::set_Domain(function<float(float)> f_top,function<float(float)> f_bottom){
    printf("wow\n");
//    dom.SetZero();
    cout << size_x << " and " << size_y << endl;
    cout << "dom size = " << dom.size() << endl;
    for(int k=1; k<=size_x*size_y-1; k++) {
        cout << "k = " << k << endl;
        dom[k] = false;
        if( ( k/size_x > f_top(k%size_x) ) && ( k/size_x < size_y - f_bottom(k%size_x) ) ){
            cout << "k true = " << k << endl;
            dom[k] = true;
        }
    }                 
}

/* codigo em python dado pelo Cosme
for i in range(200):
    for j in range(200):
#        if j>0.2*i and j< 200-0.3*i  : D[i,j]=1
#        if j>5+15*(np.tanh(i-100)+1) and j< 180-20*np.cos(i/20.)  : D[i,j]=1
        if j > 5 + (0 if i>100 else 30) and j < 180 - 20 * np.cos(i / 20.): D[i, j] = 1
*/
