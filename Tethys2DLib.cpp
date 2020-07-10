// 2D version

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>

#include <iomanip>   


#include "TethysLib.h"
#include "Tethys2DLib.h"
#include <H5Cpp.h>

using namespace H5;
using namespace std;


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

Fluid2D::Fluid2D(int sizeNx, int sizeNy, float VELSND, float VISCO) : TETHYSBase{sizeNx,sizeNy,2}{		
	Nx = sizeNx;
	Ny = sizeNy;
	vel_snd =VELSND;
	kin_vis =VISCO;
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
	// main grid variables Nx*Ny
	den 		= new float[Nx*Ny](); 
	velX 		= new float[Nx*Ny](); 
	velY 		= new float[Nx*Ny](); 
	flxX 		= new float[Nx*Ny](); 
	flxY 		= new float[Nx*Ny](); 
	curX 		= new float[Nx*Ny](); 
	curY 		= new float[Nx*Ny]();
	vel_snd_arr	= new float[Nx*Ny]();

    lap_flxX    = new float[Nx*Ny](); //new grids for the laplacians
    lap_flxY    = new float[Nx*Ny](); //in fact they could be smaller but thiw way they are just 0 at the borders who do not evolve

	// 1st Aux. Grid variables (Nx-1)*(Ny-1)
	den_mid		= new float[(Nx-1)*(Ny-1)]();  
	flxX_mid	= new float[(Nx-1)*(Ny-1)]();
	flxY_mid	= new float[(Nx-1)*(Ny-1)]();

}
	
Fluid2D::~Fluid2D(){
	delete [] den;
	delete [] velX;
	delete [] velY;
	delete [] flxX;
	delete [] flxY;
	delete [] curX;
	delete [] curY;
	delete [] den_mid;
	delete [] flxX_mid;
	delete [] flxY_mid;
    delete [] lap_flxX;
    delete [] lap_flxY;
    delete [] vel_snd_arr;
}



void Fluid2D::SetSound(){ 
	for(int i = 0; i<Nx  ;i++){
		for(int j=0; j<Ny ; j++){
			//vel_snd_arr[i+j*Nx]=SoundVelocityAnisotropy(i,dx,j,dy, vel_snd);
			vel_snd_arr[i+j*Nx]=vel_snd;
		}
	}
}


float Fluid2D::GetVelSnd(){ return vel_snd; }
void Fluid2D::SetVelSnd(float x){ vel_snd=x; }
float Fluid2D::GetKinVis(){ return kin_vis; }
void Fluid2D::SetKinVis(float x){ kin_vis=x;}
float Fluid2D::GetDx(){return dx;}
void Fluid2D::SetDx(float x){ dx=x;} 
float Fluid2D::GetDy(){return dy;}
void Fluid2D::SetDy(float x){ dy=x;}
float Fluid2D::GetDt(){return dt;}
void Fluid2D::SetDt(float x){ dt=x;}

void Fluid2D::InitialCondRand(){
  	srand (static_cast<unsigned int>(time(NULL)));
  	for (int i = 0; i < Nx; i++ ){
  		for (int j=0; j<Ny; j++){
  			float noise = (float) rand()/ (float) RAND_MAX ;
			den[i+j*Nx] = 1.0f + 0.005f*(noise-0.5f);
  		}
  	}	
}

void Fluid2D::InitialCondTEST(){
  	for (int i = 0; i < Nx; i++ ){
  		for (int j=0; j<Ny; j++){
			float densi;
			if(i>=80&&i<=120&&j>=80&&j<=120){
			densi=0.2f;
			}
			else{
			densi=0.0f;
			}
			den[i+j*Nx] = 1.0f + densi;
			velX[i+j*Nx] = 0.1f;
  		}
  	}	
}


void Fluid2D::MassFluxToVelocity(){
	for(int C=0;C<=Nx*Ny-1;C++){
		velX[C]=flxX[C]/den[C];
		velY[C]=flxY[C]/den[C];
		curX[C] =velX[C]*den[C];
		curY[C] =velY[C]*den[C];			
	}
}

//void Fluid2D::DimensionalSplittingMethod(){
//// x-sweeps runing all lines from 0<y<W i.e. excluding the boundaries
//	for(int j=1;j<Ny-1;j++){
//		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
//		for ( int i = 0; i < Nx - 1; i++ )
//		{
//			den_mid[i+j*Nx] = 0.5*( den[i+j*Nx] + den[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( DensityFluxX(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - DensityFluxX(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//			flxX_mid[i+j*Nx] = 0.5*( flxX[i+j*Nx] + flxX[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( MassFluxXFluxX(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - MassFluxXFluxX(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//			flxY_mid[i+j*Nx] = 0.5*( flxY[i+j*Nx] + flxY[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( MassFluxYFluxX(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - MassFluxYFluxX(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//		}
//		// Remaining step
//		for ( int i = 1; i < Nx - 1; i++ )
//		{
//			den[i+j*Nx] =  den[i+j*Nx]  - (dt/dx) * ( DensityFluxX(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - DensityFluxX(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//			flxX[i+j*Nx] = flxX[i+j*Nx] - (dt/dx) * ( MassFluxXFluxX(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - MassFluxXFluxX(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//			flxY[i+j*Nx] = flxY[i+j*Nx] - (dt/dx) * ( MassFluxYFluxX(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - MassFluxYFluxX(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//		}
//	}
//// y-sweeps
//	for(int i=1;i<Nx-1;i++){
//		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
//		for ( int j = 0; j < Ny - 1; j++ )
//		{
//			den_mid[i+j*Nx] = 0.5*( den[i+j*Nx] + den[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( DensityFluxY(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - DensityFluxY(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//			flxX_mid[i+j*Nx] = 0.5*( flxX[i+j*Nx] + flxX[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( MassFluxXFluxY(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - MassFluxXFluxY(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//			flxY_mid[i+j*Nx] = 0.5*( flxY[i+j*Nx] + flxY[i+1+j*Nx] )
//				- ( 0.5*dt/dx ) * ( MassFluxYFluxY(den[i+1+j*Nx],flxX[i+1+j*Nx],flxY[i+1+j*Nx],vel_snd) - MassFluxYFluxY(den[i+j*Nx],flxX[i+j*Nx],flxY[i+j*Nx],vel_snd) );
//		}
//		// Remaining step
//		for ( int j = 1; j < Ny - 1; j++ )
//		{
//			den[i+j*Nx] =  den[i+j*Nx]  - (dt/dx) * ( DensityFluxY(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - DensityFluxY(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//			flxX[i+j*Nx] = flxX[i+j*Nx] - (dt/dx) * ( MassFluxXFluxY(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - MassFluxXFluxY(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//			flxY[i+j*Nx] = flxY[i+j*Nx] - (dt/dx) * ( MassFluxYFluxY(den_mid[i+j*Nx],flxX_mid[i+j*Nx],flxY_mid[i+j*Nx],vel_snd) - MassFluxYFluxY(den_mid[i-1+j*Nx],flxX_mid[i-1+j*Nx],flxY_mid[i-1+j*Nx],vel_snd) );
//		}
//	}
//}


void Fluid2D::Richtmyer(){
    // TODO throw exception to cath NAN or INF values
		int NE,NW,SE,SW;
		float n_N, n_S ,n_E ,n_W, px_N, px_S, px_E, px_W, py_N, py_S, py_E, py_W,m_E,m_W,m_N,m_S;
		//k=i+j*Nx
		for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
            div_t divresult;
            divresult = div (ks,Nx-1);
            int j=divresult.quot;
            int i=divresult.rem;

				NE=i+1+(j+1)*Nx; //mal  ->i+1,j+1 Prin  kPrin = i+j*Nx
				NW=i+(j+1)*Nx;   //mal  ->i,j+1   Prin
				SE=i+1+j*Nx;    //mal  ->i+1,j   Prin
				SW=i+j*Nx;      //mal  ->i,j     Prin
		
				n_N = 0.5f*(den[NE]+den[NW]);
				n_S = 0.5f*(den[SE]+den[SW]);
				n_E = 0.5f*(den[NE]+den[SE]);
				n_W = 0.5f*(den[NW]+den[SW]);

				px_N = 0.5f*(flxX[NE]+flxX[NW]);
				px_S = 0.5f*(flxX[SE]+flxX[SW]);
				px_E = 0.5f*(flxX[NE]+flxX[SE]);
				px_W = 0.5f*(flxX[NW]+flxX[SW]);
				
				py_N = 0.5f*(flxY[NE]+flxY[NW]);
				py_S = 0.5f*(flxY[SE]+flxY[SW]);
				py_E = 0.5f*(flxY[NE]+flxY[SE]);
				py_W = 0.5f*(flxY[NW]+flxY[SW]);

				//posso definir aqui a "massa" nos 4 ponto s cardeais
				 m_E=pow(n_E,1.5f); // e assim sucessivamente m_W m_N m_S que depois sao reutilizadeas nos 12 fluxos
				 m_W=pow(n_W,1.5f);
                 m_N=pow(n_N,1.5f);
                 m_S=pow(n_S,1.5f);
				den_mid[ks] = 0.25f*(den[SW] + den[SE] + den[NW] + den[NE]) // How shall we include vel_snd_arr ?
								-0.5f*(dt/dx)*(
									DensityFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									DensityFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-0.5f*(dt/dy)*(
									DensityFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									DensityFluxY(n_S, px_S, py_S,m_S,vel_snd));
				flxX_mid[ks] = 0.25f*(flxX[SW] + flxX[SE] + flxX[NW] + flxX[NE])
								-0.5f*(dt/dx)*(
									MassFluxXFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									MassFluxXFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-0.5f*(dt/dy)*(
									MassFluxXFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									MassFluxXFluxY(n_S, px_S, py_S,m_S,vel_snd));
				flxY_mid[ks] = 0.25f*(flxY[SW] + flxY[SE] + flxY[NW] + flxY[NE])
								-0.5f*(dt/dx)*(
									MassFluxYFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									MassFluxYFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-0.5f*(dt/dy)*(
									MassFluxYFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									MassFluxYFluxY(n_S, px_S, py_S,m_S,vel_snd));
	//		}						
		}
		for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
            div_t divresult;
            divresult = div (kp,Nx);
            int j=divresult.quot;
            int i=divresult.rem;


			if( kp%Nx!=Nx-1 && kp%Nx!=0){
				
				NE=i+j*(Nx-1);     //mal  ->i,j Sec  kSec = i+j*(Nx-1)
				NW=i-1+j*(Nx-1);   //mal  ->i-1,j Sec
				SE=i+(j-1)*(Nx-1);  //mal  ->i,j-1 Sec
				SW=i-1+(j-1)*(Nx-1); //mal ->i-1,j-1 Sec
				
				n_N = 0.5f*(den_mid[NE]+den_mid[NW]);
				n_S = 0.5f*(den_mid[SE]+den_mid[SW]);
				n_E = 0.5f*(den_mid[NE]+den_mid[SE]);
				n_W = 0.5f*(den_mid[NW]+den_mid[SW]);

				px_N = 0.5f*(flxX_mid[NE]+flxX_mid[NW]);
				px_S = 0.5f*(flxX_mid[SE]+flxX_mid[SW]);
				px_E = 0.5f*(flxX_mid[NE]+flxX_mid[SE]);
				px_W = 0.5f*(flxX_mid[NW]+flxX_mid[SW]);
				
				py_N = 0.5f*(flxY_mid[NE]+flxY_mid[NW]);
				py_S = 0.5f*(flxY_mid[SE]+flxY_mid[SW]);
				py_E = 0.5f*(flxY_mid[NE]+flxY_mid[SE]);
				py_W = 0.5f*(flxY_mid[NW]+flxY_mid[SW]);

                //posso definir aqui a "massa" nos 4 ponto s cardeais
                 m_E=pow(n_E,1.5f); // e assim sucessivamente m_W m_N m_S que depois sao reutilizadeas nos 12 fluxos
                 m_W=pow(n_W,1.5f);
                 m_N=pow(n_N,1.5f);
                 m_S=pow(n_S,1.5f);

				den[kp] = den[kp]
								-(dt/dx)*(
									DensityFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									DensityFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-(dt/dy)*(
									DensityFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									DensityFluxY(n_S, px_S, py_S,m_S,vel_snd));
				flxX[kp] = flxX[kp]
								-(dt/dx)*(
									MassFluxXFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									MassFluxXFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-(dt/dy)*(
									MassFluxXFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									MassFluxXFluxY(n_S, px_S, py_S,m_S,vel_snd));
				flxY[kp] = flxY[kp]
								-(dt/dx)*(
									MassFluxYFluxX(n_E, px_E, py_E,m_E,vel_snd)-
									MassFluxYFluxX(n_W, px_W, py_W,m_W,vel_snd))
								-(dt/dy)*(
									MassFluxYFluxY(n_N, px_N, py_N,m_N,vel_snd)-
									MassFluxYFluxY(n_S, px_S, py_S,m_S,vel_snd));
			}
		}	
	
}
		


void Fluid2D::CFLCondition(){ 
		dx = lengX / ( float ) ( Nx - 1 );
		dy = lengY / ( float ) ( Ny - 1 );
		dt = dx/10.0f;
}



float  Fluid2D::DensityFluxX(float n,float flxX, float flxY,float mass, float S){
	float f1;
	f1 = flxX;
	return f1;		
}
float  Fluid2D::DensityFluxY(float n,float flxX, float flxY,float mass, float S){
	float f1;
	f1 = flxY;
	return f1;		
}
float  Fluid2D::DensitySource(float n,float velX, float velY, float S){
	float Q1 =0.0f;
	return Q1;
}
float  Fluid2D::MassFluxXFluxX(float n,float flxX, float flxY,float mass, float S){
	float f2;
	f2 = flxX*flxX/n +n; 
	return f2;
}
float  Fluid2D::MassFluxXFluxY(float n,float flxX, float flxY,float mass, float S){
	float f2; 
	f2 = flxX*flxY/n;
	return f2;
}
float  Fluid2D::MassFluxYFluxX(float n,float flxX, float flxY,float mass, float S){
	float f3;
	f3 = flxX*flxY/n;
	return f3;
}
float  Fluid2D::MassFluxYFluxY(float n,float flxX, float flxY,float mass, float S){
	float f3;
	f3 = flxY*flxY/n + n;
	return f3;
}
float  Fluid2D::MassFluxXSource(float n,float flxX, float flxY, float S){
	float Q2 =0.0f;
	return Q2;
}
float  Fluid2D::MassFluxYSource(float n,float flxX, float flxY, float S){
	float Q3 =0.0f;
	return Q3;
}	

void Fluid2D::SetFileName(){
	char buffer [50];
	sprintf (buffer, "S=%.2fvis=%.2f", vel_snd, kin_vis);
	file_infix = buffer;
}

void Fluid2D::CreateFluidFile(){
	this->SetFileName();
	std::string previewfile = "preview_2D_" + file_infix + ".dat" ;
	data_preview.open (previewfile);
	data_preview << scientific; 
}

void Fluid2D::WriteFluidFile(float t){
int j=Ny/2;
data_preview <<t<<"\t"<< den[Nx-1+j*Nx] <<"\t"<< flxX[Nx-1+j*Nx] <<"\t"<< den[0+j*Nx] <<"\t" << flxX[0+j*Nx] <<"\n";
}



void Fluid2D::SetSimulationTime(){
	Tmax=5.0f+0.02f*vel_snd+20.0f/vel_snd;
}

GrapheneFluid2D::GrapheneFluid2D(int sizeNx,int sizeNy,float VELSND, float FERMI,float VISCO,float COL): Fluid2D(sizeNx,sizeNy, VELSND, VISCO){
	vel_fer =FERMI;							
	col_freq =COL; 
	char buffer [50];
	sprintf (buffer, "S=%.2fvF=%.2fvis=%.2fl=%.2f", vel_snd, vel_fer, kin_vis, col_freq);
	file_infix = buffer;
}


void GrapheneFluid2D::SetSimulationTime(){
	float s;
	s=this->GetVelSnd();
	this->SetTmax(5.0f+0.02f*s+20.0f/s);
}

void GrapheneFluid2D::MassFluxToVelocity(){
	for(int C=0;C<=Nx*Ny-1;C++){
		velX[C]=flxX[C]*pow(den[C],-1.5f);
		velY[C]=flxY[C]*pow(den[C],-1.5f);
		curX[C] = velX[C]*den[C];
		curY[C] = velY[C]*den[C];			
	}
}


void GrapheneFluid2D::SetVelFer(float x){ vel_fer=x;	}
float GrapheneFluid2D::GetVelFer(){ return vel_fer;  }
void GrapheneFluid2D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid2D::GetColFreq(){ return col_freq; }


void GrapheneFluid2D::CFLCondition(){ // Eventual redefinition 
	dx = lengX / ( float ) ( Nx - 1 );
	dy = lengY / ( float ) ( Ny - 1 );					
	//dt = 2.4/(vel_snd*sqrt(25.0/(dx*dx)+16.0/(dy*dy)));
	float lambda;
	if(vel_snd<0.36f*vel_fer){
		lambda=1.2f*vel_fer;
	}else{
		lambda=1.97f*vel_snd + 0.5f*vel_fer;
	}
	dt = dx/lambda;				
}	

float  GrapheneFluid2D::DensityFluxX(float n,float flxX, float flxY,float mass, float S){ // Double Please Review this João
	float f1;
	f1 = flxX/sqrt(n);
	return f1;		
}
float  GrapheneFluid2D::DensityFluxY(float n,float flxX, float flxY,float mass, float S){
	float f1;
	f1 = flxY/sqrt(n);
	return f1;		
}
// 27% of cpu usage of them 22.6% are int the pow function
float  GrapheneFluid2D::MassFluxXFluxX(float n,float flxX, float flxY,float mass, float S){
	float f2;
	f2 = flxX*flxX/mass +vel_fer*vel_fer*mass/3.0f+0.5f*S*S*n*n;
	return f2;
}
float  GrapheneFluid2D::MassFluxXFluxY(float n,float flxX, float flxY,float mass, float S){
	float f2; 
	f2 = flxX*flxY/mass;
	return f2;
}
float  GrapheneFluid2D::MassFluxYFluxX(float n,float flxX, float flxY,float mass, float S){
	float f3;
	f3 = flxX*flxY/mass;
	return f3;
}
float  GrapheneFluid2D::MassFluxYFluxY(float n,float flxX, float flxY,float mass, float S){
	float f3;
	f3 = flxY*flxY/mass + vel_fer*vel_fer*mass/3.0f+0.5f*S*S*n*n;
	return f3;
}




// Pedro: para ja nao vamos incluir sources 
float  GrapheneFluid2D::DensitySource(float n,float flxX, float flxY, float S){
	float Q1 =0;
	return Q1;
}
float  GrapheneFluid2D::MassFluxXSource(float n,float flxX, float flxY, float S){
	float Q2 =0;
	return Q2;
}
float  GrapheneFluid2D::MassFluxYSource(float n,float flxX, float flxY, float S){
	float Q3 =0;
	return Q3;
}	


void GrapheneFluid2D::MagneticSource(){
	float px0,py0,sqrtn0;
	float Wc=10.0;
    for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
        if( kp%Nx!=Nx-1 && kp%Nx!=0){
			sqrtn0=sqrt(den[kp]);
			px0=flxX[kp];
			py0=flxY[kp];
			flxX[kp]=px0*cos(Wc*dt/sqrtn0)-py0*sin(Wc*dt/sqrtn0);
			flxY[kp]=px0*sin(Wc*dt/sqrtn0)+py0*cos(Wc*dt/sqrtn0);
		}
	}		
}

void GrapheneFluid2D::ViscosityFTCS() {
    int N, S, E, W;
    float mass_den_C, mass_den_N, mass_den_S, mass_den_E, mass_den_W;
    //calculate laplacians
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
        div_t divresult;
        divresult = div(kp, Nx);
        int j = divresult.quot;
        int i = divresult.rem;
        if (kp % Nx != Nx - 1 && kp % Nx != 0) {
            N = i + (j + 1) * Nx;
            S = i + (j - 1) * Nx;
            E = i + 1 + j * Nx;
            W = i - 1 + j * Nx;
            mass_den_C = pow(den[kp], 1.5f);
            mass_den_N = pow(den[N], 1.5f);
            mass_den_S = pow(den[S], 1.5f);
            mass_den_E = pow(den[E], 1.5f);
            mass_den_W = pow(den[W], 1.5f);
            lap_flxX[kp] =
                    (-4.0 * flxX[kp] / mass_den_C + flxX[N] / mass_den_N + flxX[S] / mass_den_S + flxX[E] / mass_den_E +
                     flxX[W] / mass_den_W) / (dx * dx);
            lap_flxY[kp] =
                    (-4.0 * flxY[kp] / mass_den_C + flxY[N] / mass_den_N + flxY[S] / mass_den_S + flxY[E] / mass_den_E +
                     flxY[W] / mass_den_W) / (dx * dx);
        }
    }

    //FTCS algorithm
    float old_px,old_py;
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) { //correr a grelha principal evitando as fronteiras
        if (kp % Nx != Nx - 1 && kp % Nx != 0) {
            old_px=flxX[kp];
            old_py=flxY[kp];
            flxX[kp] = old_px + dt*( kin_vis*lap_flxX[kp] );
            flxY[kp] = old_py + dt*( kin_vis*lap_flxY[kp] );
        }
    }
}
void GrapheneFluid2D::SourceFTCS(){
    float px0,py0,sqrtn0;
    float Wc=10.0;
    for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
        if( kp%Nx!=Nx-1 && kp%Nx!=0){
            sqrtn0=sqrt(den[kp]);
            px0=flxX[kp];
            py0=flxY[kp];
            flxX[kp]=  px0 - 0.5f*dt*Wc*py0/sqrtn0;
            flxY[kp]=  py0 + 0.5f*dt*Wc*px0/sqrtn0;
        }
    }
}

void GrapheneFluid2D::WriteAtributes(){
	const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
	const IntType        hdf5_int(PredType::NATIVE_INT);
	int total_steps= static_cast<int>(Tmax / dt);
	//Create the data space for the attribute.
	hsize_t dim_atr[1] = { 1 };
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = grp_dat->createAttribute( "S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = grp_dat->createAttribute( "Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_kin_vis = grp_dat->createAttribute( "Kinetic viscosity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = grp_dat->createAttribute( "Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = grp_dat->createAttribute( "Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = grp_dat->createAttribute( "Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = grp_dat->createAttribute( "Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_space_points = grp_dat->createAttribute( "Number of spatial points", hdf5_int, atr_dataspace);
	Attribute atr_num_time_steps = grp_dat->createAttribute( "Number of time steps", hdf5_int, atr_dataspace);
	// Write the attribute data.
	atr_vel_snd.write( hdf5_float, &vel_snd);
	atr_vel_fer.write( hdf5_float, &vel_fer);
	atr_col_freq.write(hdf5_float, &col_freq);
	atr_kin_vis.write(hdf5_float, &kin_vis); 
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_num_space_points.write( hdf5_int, &Nx);
	atr_total_time.write( hdf5_float, &Tmax);
	atr_num_time_steps.write(hdf5_int, &total_steps);
	// Close the attributes.
	atr_num_time_steps.close();
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_kin_vis.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
}
