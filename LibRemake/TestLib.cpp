
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
#include "TestLib.h"





Fluid1D::Fluid1D(int sizeN){		
	Nx = sizeN;
	den     = new float[sizeN]();
	vel     = new float[sizeN]();
	cur     = new float[sizeN]();
	den_cor = new float[sizeN]();
	vel_cor = new float[sizeN]();
	cur_cor = new float[sizeN]();
	den_mid = new float[sizeN-1]();
	vel_mid = new float[sizeN-1]();
	vel_snd_arr = new float[sizeN-1]();
}	
	
Fluid1D::~Fluid1D(){			
	delete [] den;
	delete [] vel ;
	delete [] cur ;
	delete [] den_mid ;
	delete [] vel_mid ;
	delete [] den_cor ;
	delete [] vel_cor ;
	delete [] cur_cor ;
	delete [] vel_snd_arr ;
}

float  Fluid1D::DensityFlux(float n,float v,float S){
	float f1;
	f1 = n*v;
	return f1;		
}
float  Fluid1D::VelocityFlux(float n,float v,float S){
	float f2;
	f2 = 0.5*v*v + n; 
	return f2;
}
float  Fluid1D::DensitySource(float n,float v,float S){
	return 0;
}
float  Fluid1D::VelocitySource(float n,float v,float S){
	return 0;
}	

void Fluid1D::CFLCondition(){
		dx = leng / ( float ) ( Nx - 1 );
		dt = dx/10.0;
}
		
void Fluid1D::SetSound(){
	for(int i = 0; i<Nx  ;i++){
		vel_snd_arr[i]=SoundVelocityAnisotropy(i,dx,vel_snd);
	}
}
		
void Fluid1D::InitialCondRand(){
  srand (time(NULL));   
  for (int i = 0; i < Nx; i++ )
  {
		float noise = (float) rand()/ (float) RAND_MAX ;
		den[i] = 1.0 + 0.005*(noise-0.5);
  }	
}


void Fluid1D::SetVelSnd(float x){ vel_snd=x; }
float Fluid1D::GetVelSnd(){ return vel_snd; }
float Fluid1D::GetDx(){return dx;}
void Fluid1D::SetDx(float x){ dx=x;}
float Fluid1D::GetDt(){return dt;}
void Fluid1D::SetDt(float x){ dt=x;}
int Fluid1D::SizeX(){ return Nx; }

void Fluid1D::Smooth(int width){
	AverageFilter( den ,den_cor, Nx, width);	
	AverageFilter( vel ,vel_cor, Nx, width);
	AverageFilter( cur ,cur_cor, Nx , width);
}

void Fluid1D::Richtmyer(){
    	//
		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i < Nx - 1; i++ )
		{
			den_mid[i] = 0.5*( den[i] + den[i+1] )
				- ( 0.5*dt/dx ) * ( DensityFlux(den[i+1],vel[i+1],vel_snd_arr[i]) - DensityFlux(den[i],vel[i],vel_snd_arr[i]) ) 
				+ ( 0.5*dt    ) * DensitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),vel_snd_arr[i]) ;
			vel_mid[i] = 0.5*( vel[i] + vel[i+1] )
				- ( 0.5*dt/dx ) * ( VelocityFlux(den[i+1],vel[i+1],vel_snd_arr[i]) - VelocityFlux(den[i],vel[i],vel_snd_arr[i]) ) 
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),vel_snd_arr[i]) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],vel_snd_arr[i]) - DensityFlux(den_mid[i-1],vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * DensitySource(den[i],vel[i],vel_snd_arr[i]);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],vel_snd_arr[i]) - VelocityFlux(den_mid[i-1],vel_mid[i-1],vel_snd_arr[i]) )
							+  dt * VelocitySource(den[i],vel[i],vel_snd_arr[i]);
			cur[i] = vel[i]*den[i];
		}
} 




		float GrapheneFluid1D::DensityFlux(float n,float v,float S){
			float f1;
			f1 = n*v;
			return f1;			
		}
		float GrapheneFluid1D::VelocityFlux(float n,float v,float S){
			float f2;
			f2 = 0.25*v*v + vel_fer*vel_fer*0.5*log(n) + 2*S*S*sqrt(n); 
			return f2;			
		}
		float GrapheneFluid1D::DensitySource(float n,float v,float S){
			float Q1=0.0;
			return Q1;				
		}
		float GrapheneFluid1D::VelocitySource(float n,float v,float S){
			float Q2=0.0;
			Q2=-1.0*col_freq*(v-1);
			return Q2;			
		}
		



void GrapheneFluid1D::SetVelFer(float x){ vel_fer=x;	}
float GrapheneFluid1D::GetVelFer(){ return vel_fer;  }
void GrapheneFluid1D::SetColFreq(float x){ col_freq=x; }
float GrapheneFluid1D::GetColFreq(){ return col_freq; }

	 
void GrapheneFluid1D::CFLCondition(){
	dx = leng / ( float ) ( Nx - 1 );
					
	if(vel_fer<10 && (vel_snd-vel_fer) <= 3)
		dt = 0.5 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<10 && (vel_snd-vel_fer<= 10 - vel_fer))
		dt = 1.5 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<15 && (vel_snd-vel_fer<= 5))
		dt = 2 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else if (vel_fer<30 && (vel_snd-vel_fer<= 3))
		dt = 3 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
	else
		dt = 4 * dx / (2*vel_snd+sqrt(3*vel_fer*vel_fer + 24*vel_snd*vel_snd));
}	


void GrapheneFluid1D::BoundaryCond(int type){
	
		/*---------------*\
	| Free        | 1 |
	| Periodic 	  | 2 |
	| DS Boundary | 3 | 
	| DS+Driving  | 4 | 
	\*---------------*/
	switch(type){
		case 1 : den[0] = den[1];
				 den[Nx-1] = den[Nx-2];
				 vel[0] = vel[1];
				 vel[Nx-1] = vel[Nx-2];		
			break;
		case 2 : den[0] = den[Nx-2];
				 den[Nx-1] = den[1];
				 vel[0] = vel[Nx-2];
				 vel[Nx-1] = vel[1];	
			break;
		case 3 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[Nx-1] = den[Nx-2];
				 vel[Nx-1] = 1.0/den[Nx-1];			
			break;	
		case 4 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[Nx-1] = den[Nx-2];
				 vel[Nx-1] = (1.0 + 0.75*vel[1]*den[1])/den[Nx-1];			
			break;		
		default : den[0] = 1.0;
				  vel[0] = vel[1];
				  den[Nx-1] = den[Nx-2];
				  vel[Nx-1] = 1.0/den[Nx-1];			
	}

	
}

float SoundVelocityAnisotropy(float i, float dx,float S){
	return S;
}


void AverageFilter(float * vec_in, float * vec_out, int size , int width ){
			for ( int i = 0; i < size; i++ ){		
				if(i>=width &&i<=size-1-width){
					for(int k = i-width; k <= i+width;k++){			
						vec_out[i] += vec_in[k]; 
					}	
					vec_out[i] = vec_out[i]/(2.0*width+1.0);
					}
				else{
					vec_out[i] =	vec_in[i] ;
				}	
			}	
		}

