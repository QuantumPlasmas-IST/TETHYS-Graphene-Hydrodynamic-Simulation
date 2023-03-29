#include "includes/DomainLib.h"

void Domain::fill_Domain(Fluid2D &fluid_class){
    int Nx = fluid_class.SizeX();
    int Ny = fluid_class.SizeY();
    for(int i=1; i<=Nx*Ny-1; i++) {
        if( ( i%Nx > Domain::bottom_margin(i/Nx) ) && ( i%Nx < Ny - Domain::top_margin(i/Nx) ) ){
            dom[i/Nx , i%Nx] = true;
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
