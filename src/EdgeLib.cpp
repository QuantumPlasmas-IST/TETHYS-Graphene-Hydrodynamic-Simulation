#include "EdgeLib.h"
#include "DomainLib.h"

Edge::Edge(){}

Edge::Edge(int Nx, int Ny /*,bool * domn*/){
    size_x = Nx;
    size_y = Ny;
//    D.dom = domn;
    D.set_size_x(Nx);
    D.set_size_y(Ny);
    edg.resize(Nx*Ny);
//    edgint.resize(Ny*Nx);
}

Edge::~Edge() = default;

void Edge::set_size_x(int Nx){
    size_x = Nx;
    D.set_size_x(Nx);
}

void Edge::set_size_y(int Ny){
    size_y = Ny;
    D.set_size_y(Ny);
}

void Edge::set_Edge(){
    int edgctr = 0;
    for(int k=1; k<=size_x*size_y-1; k++) {
//        cout << "k = " << k << endl;
//        printf("i/size_x =%d",k/size_x);
//        cout << "D.dom[" << k << "] = " << D.dom[k] << endl;
        if(D.dom[k] == 0){
            for(int di = -1; di <= 1; di++){
                for(int dj = -1; dj <= 1; dj++){
//                    cout << "D.dom[" << k + di + size_x*dj << "] = " << D.dom[k + di + size_x*dj] << endl;
                    if( D.dom[k + di + size_x*dj] == 1){
//                        cout << "edgctr =" << edgctr << endl;
//                        cout << "edgint =" << edgint[edgctr] << endl;
                        edg[k] = 1;
                        edgctr ++;
                        edgint.push_back(k);
                    //    edgint[edgctr]=k;
                    //    edgctr++; 
                    }
                }
            }
        }
    }
    cout << "edgctr = " << edgctr << endl;                 
}

/* codigo em python dado pelo Cosme
for i in range(200):
    for j in range(200):
#        if j>0.2*i and j< 200-0.3*i  : D[i,j]=1
#        if j>5+15*(np.tanh(i-100)+1) and j< 180-20*np.cos(i/20.)  : D[i,j]=1
        if j > 5 + (0 if i>100 else 30) and j < 180 - 20 * np.cos(i / 20.): D[i, j] = 1


for i in range(1,199,1):
    for j in range(1,199,1):
        if D[i,j]==0 :
            for di in {-1,0,1}:
                for dj in {-1, 0, 1}:
                    if D[i+di, j+dj] == 1: E[i,j]=1*/


