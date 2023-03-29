#include "EdgeLib.h"

void Edge::condition_Edge(Fluid2D &fluid_class){
    int Nx = fluid_class.SizeX();
    int Ny = fluid_class.SizeY();
    
    for(int i=1; i<=Nx*Ny-1; i++) {
        printf("i/Nx =%d",i/Nx);
        if(dom[i/Nx,i%Nx] == 0){
            for(int di = -1; di <= 1; di++){
                for(int dj = -1; dj <= 1; dj++){
                    if( dom[i/Nx + di,i%Nx + di] == 0){
                        if( dom[i/Nx + dj,i%Nx + dj] == 1){
                            edg[i/Nx,i%Nx] == 1;
                        }
                    }
                }
            }
        }
    }                 
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


