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
#include "EM_fields.h"


#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

#ifndef C_SPEED
#    define C_SPEED 1000.0
#endif


using namespace std;

float RetardedTime(float time, float x , float y , float z, float X , float Y , float Z ){
	float tr;
	tr = time - CartDistance( x, y, z, X, Y, Z)/C_SPEED;
	if(tr>=0){
		return tr;		
	}
	else{
		return time;	
	}
}


void JefimenkoEMField(int XDIM, int YDIM, float dx, float dy, float dt, float Xpos, float Ypos, float Zpos,  float ** rho, float ** rho_dot, float ** cur, float ** cur_dot, float Time , float  * E_out , float  * B_out, float  * S_out   ){
	//float q =-1.0;
	int k_retard;
	float R_norm;
	
	float x0, y0;
	x0=0.0;
	y0=-0.5;
	float x,y,z;
	z=0.0;
	
	int N,M;
	N = XDIM;
	M = YDIM;
	
	float SumX_e0=0.0; 
	float SumY_e0=0.0; 
	float SumDiag_e0=0.0;
	float Corners_e0=0.0;
			
	float SumX_cur=0.0; 
	float SumY_cur=0.0; 
	float SumDiag_cur=0.0;
	float Corners_cur=0.0;
			
	float SumX_by=0.0; 
	float SumY_by=0.0; 
	float SumDiag_by=0.0;
	float Corners_by=0.0;
	
	float E0 = 0.0;
	float Cur= 0.0;
	float B0 = 0.0;	
	
	for(int i=1;i<=N-1;i++){
		/*y = 0 */
		x = x0 + i*dx;
		y = y0 + 0*dy;;
		k_retard = nearbyint( RetardedTime(Time , x , y , z, Xpos, Ypos , Zpos )/dt);	
		R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );		
		SumX_e0  +=        ( rho[ k_retard][i]/pow(R_norm,3) + rho_dot[ k_retard][i]/(pow(R_norm,2)*C_SPEED) );
		SumX_cur += -1.0*x*( rho[ k_retard][i]/pow(R_norm,3) + rho_dot[ k_retard][i]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][i]/(R_norm*C_SPEED*C_SPEED);
		SumX_by  +=        ( cur[ k_retard][i]/pow(R_norm,3) + cur_dot[ k_retard][i]/(pow(R_norm,2)*C_SPEED) );
		/*y = M */
		x = x0 + i*dx;
		y = y0 + M*dy;
		k_retard = nearbyint( RetardedTime(Time , x , y , z, Xpos, Ypos , Zpos )/dt);	
		R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );	
		SumX_e0  +=        ( rho[k_retard][i]/pow(R_norm,3) + rho_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED) );
		SumX_cur += -1.0*x*( rho[k_retard][i]/pow(R_norm,3) + rho_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][i]/(R_norm*C_SPEED*C_SPEED);
		SumX_by  +=        ( cur[k_retard][i]/pow(R_norm,3) + cur_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED) );
	}		
				
	for(int j=1;j<=M-1;j++){
		/*x = 0*/
		x = x0 + 0*dx;
		y = y0 + j*dy;
		k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
		R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
		SumY_e0  +=        ( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) );
		SumY_cur += -1.0*x*( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][0]/(R_norm*C_SPEED*C_SPEED);
		SumY_by  +=        ( cur[k_retard][0]/pow(R_norm,3) + cur_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) ); 
		/*x = N*/
		x = x0 + N*dx;
		y = y0 + j*dy;
		k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
		R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
		SumY_e0  +=        ( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) );
		SumY_cur += -1.0*x*( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][N]/(R_norm*C_SPEED*C_SPEED);
		SumY_by  +=        ( cur[k_retard][N]/pow(R_norm,3) + cur_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) ); 
	}		
			
	for(int j=1;j<=M-1;j++){
		for(int i=1;i<=N-1;i++){
			x = x0 + i*dx;
			y = y0 + j*dy;
			k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
			R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
			SumDiag_e0  +=        ( rho[k_retard][i]/pow(R_norm,3) + rho_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED));
			SumDiag_cur += -1.0*x*( rho[k_retard][i]/pow(R_norm,3) + rho_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED)) - cur_dot[k_retard][i]/(R_norm*C_SPEED);
			SumDiag_by  +=        ( cur[k_retard][i]/pow(R_norm,3) + cur_dot[k_retard][i]/(pow(R_norm,2)*C_SPEED));
		}			
	}	
	// (i=0,j=0)
	x = x0 + 0*dx;
	y = y0 + 0*dy;
	k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
	R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
	Corners_e0  +=        ( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) );
	Corners_cur += -1.0*x*( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][0]/(R_norm*C_SPEED*C_SPEED);
	Corners_by  +=        ( cur[k_retard][0]/pow(R_norm,3) + cur_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED));
	// (i=N,j=0)
	x = x0 + N*dx;
	y = y0 + 0*dy;
	k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
	R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
	Corners_e0  +=        ( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) );
	Corners_cur += -1.0*x*( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) ) - cur_dot[k_retard][N]/(R_norm*C_SPEED*C_SPEED);
	Corners_by  +=        ( cur[k_retard][N]/pow(R_norm,3) + cur_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED));
	// (i=0,j=M)
	x = x0 + 0*dx;
	y = y0 + M*dy;
	k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
	R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
	Corners_e0  +=        ( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) );		
	Corners_cur += -1.0*x*( rho[k_retard][0]/pow(R_norm,3) + rho_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][0]/(R_norm*C_SPEED*C_SPEED);		
	Corners_by  +=        ( cur[k_retard][0]/pow(R_norm,3) + cur_dot[k_retard][0]/(pow(R_norm,2)*C_SPEED));
	// (i=N,j=M)
	x = x0 + N*dx;
	y = y0 + M*dy;
	k_retard = nearbyint( RetardedTime(Time, x , y , z, Xpos, Ypos , Zpos )/dt);	
	R_norm   =  CartDistance( x , y , z,  Xpos, Ypos , Zpos );
	Corners_e0  +=        ( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) );		 
	Corners_cur += -1.0*x*( rho[k_retard][N]/pow(R_norm,3) + rho_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED) )- cur_dot[k_retard][N]/(R_norm*C_SPEED*C_SPEED);			 
	Corners_by  +=        ( cur[k_retard][N]/pow(R_norm,3) + cur_dot[k_retard][N]/(pow(R_norm,2)*C_SPEED));
					
	E0  = 0.25*dx*dy*( Corners_e0 + 2.0*SumX_e0 + 2.0*SumY_e0 + 4.0*SumDiag_e0);		
	Cur = 0.25*dx*dy*( Corners_cur + 2.0*SumX_cur + 2.0*SumY_cur + 4.0*SumDiag_cur);					
	B0  = 0.25*dx*dy*( Corners_by + 2.0*SumX_by + 2.0*SumY_by + 4.0*SumDiag_by);		
			
	E_out[0] = E0*Xpos + Cur;
	E_out[1] = 0.0;
	E_out[2] = E0*Zpos;
	
	B_out[0] = 0.0;
	B_out[1] = -1.0*B0*Zpos;
	B_out[2] = 0.0;
		
	S_out[0] = -E_out[2]*B_out[1];
	S_out[1] = 0.0;
	S_out[2] = E_out[0]*B_out[1];
}



float TotalElectricDipole(int N,float dx, float * den){
	float dipole=0.0;
	
	for(int j=1;j<N/2;j++){	
		dipole += dx*(2*j-2)*den[2*j-2] + 4*dx*(2*j-1)*den[2*j-1] + dx*(2*j)*den[2*j];
	}
	dipole = dipole*dx/3.0;
	return dipole;
}

float TotalCurrent(int N,float dx, float * den,float * vel){	
	float cur_total=0.0;
	
	for(int j=1;j<N/2;j++){	
		cur_total += den[2*j-2]*vel[2*j-2] + 4*den[2*j-1]*vel[2*j-1] + den[2*j]*vel[2*j];
	}
	cur_total = cur_total*dx/3.0;
	return cur_total;
}
float KineticEnergy(int N,float dx, float * den, float * vel){
	float kin = 0.0;
	
	for(int j=1;j<N/2;j++){
		kin +=  0.5*vel[2*j-2]*vel[2*j-2]*den[2*j-2] + 4*0.5*vel[2*j-1]*vel[2*j-1]*den[2*j-1] + 0.5*vel[2*j]*vel[2*j]*den[2*j];
	}
	return kin*dx/3.0;
}



float DtElectricDipole(int N,float dx, float * cur){
	float dipole_deriv=0.0;
	
	for(int j=1;j<N/2;j++){
		dipole_deriv += cur[2*j-2] + 4*cur[2*j-1] + cur[2*j];
	}
	dipole_deriv = dipole_deriv*dx/3.0;
	return dipole_deriv;
}


float CartDistance(float x , float y , float z, float X , float Y , float Z ){
	float rsquared;
	rsquared = pow(X-x,2) + pow(Y-y,2) + pow(Z-z,2);
	return sqrt(rsquared);
}
