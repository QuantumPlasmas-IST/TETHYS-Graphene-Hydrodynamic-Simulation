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

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif


#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

// Time prefactor in pico seconds
#ifndef PRE_TIM
#    define PRE_TIM 3.3333333333333
#endif

// Power prefactor in pico Watts 
#ifndef PRE_POW
#    define PRE_POW 9.53395
#endif

// Kinetic Energy prefactor in atto Joules
#ifndef PRE_KIN
#    define PRE_KIN 16.8226
#endif

// Current prefactor in mili Amperes
#ifndef PRE_CUR
#    define PRE_CUR 0.480653
#endif

// VDS prefactor in Volts
#ifndef PRE_VDS
#    define PRE_VDS 18.0951
#endif

#ifndef C_SPEED
#    define C_SPEED 1000.0
#endif



using namespace std;


float RealFreq(float sound, float vel, float L, int mode){
	float real;
	if (fabs(vel) < sound ){
		mode = 2*mode-1;
	 }
	else{
		mode = 2*mode;
		}
	real =  fabs(sound*sound-vel*vel) * MAT_PI * mode / (2.0 * L * sound );
	return real;
}
	
	
float ImagFreq(float sound, float vel, float L){
	float imag;
	imag =  (sound*sound-vel*vel) * log(fabs( (sound+vel)/(sound-vel) )) / (2.0 * L * sound );
	return imag;
}	

void BoundaryCond(int type, int N, float * den, float * vel ){
	/*---------------*\
	| Free        | 1 |
	| Periodic 	  | 2 |
	| DS Boundary | 3 | 
	| DS+Driving  | 4 | 
	\*---------------*/
	switch(type){
		case 1 : den[0] = den[1];
				 den[N-1] = den[N-2];
				 vel[0] = vel[1];
				 vel[N-1] = vel[N-2];		
			break;
		case 2 : den[0] = den[N-2];
				 den[N-1] = den[1];
				 vel[0] = vel[N-2];
				 vel[N-1] = vel[1];	
			break;
		case 3 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[N-1] = den[N-2];
				 vel[N-1] = 1.0/den[N-1];			
			break;	
		case 4 : den[0] = 1.0;
				 vel[0] = vel[1];
				 den[N-1] = den[N-2];
				 vel[N-1] = (1.0 + 0.75*vel[1]*den[1])/den[N-1];			
			break;		
		default : den[0] = 1.0;
				  vel[0] = vel[1];
				  den[N-1] = den[N-2];
				  vel[N-1] = 1.0/den[N-1];			
	}
}


void InitialCondSine(int N, float dx,  float * den, float * vel){
  float L = (N-1)*dx;
  float x;
   
  for (int i = 0; i < N; i++ )
  {
		x = (float)i * dx;
		den[i] = 1.0 + 0.05*sin ( MAT_PI * x / L );
		vel[i] = 0.0;	
  }
}

void InitialCondRand(int N, float dx,  float * den, float * vel){
  srand (time(NULL)); 
  float noise; 
   
  for (int i = 0; i < N; i++ )
  {
		noise = (float) rand()/ (float) RAND_MAX ;
		den[i] = 1.0 + 0.005*(noise-0.5);
		vel[i] = 0.0;
  }
}



float DtElectricDipole(int N,float dx, float * cur){
	float dipole_deriv=0.0;
	
	for(int j=1;j<N/2;j++){
		dipole_deriv += cur[2*j-2] + 4*cur[2*j-1] + cur[2*j];
	}
	dipole_deriv = dipole_deriv*dx/3.0;
	return dipole_deriv;
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

float RootMeanSquare(int N, float dt, float * f){
	float rms=0.0;
	
	for(int j=1;j<N/2;j++){
		rms += f[2*j-2]*f[2*j-2] + 4*f[2*j-1]*f[2*j-1] + f[2*j]*f[2*j];
	}
	rms = rms*dt/3.0;
	rms = sqrt( rms/(N*dt)  );
	return rms;	
}

float SignalAverage(int N, float dt, float * f){
	float avg=0.0;
	
	for(int j=1;j<N/2;j++){
		avg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	avg = avg*dt/3.0;
	avg = avg/(N*dt);
	return avg;	
}

float GaussKernel(int position , float t){
	float g;
	
	g = exp(-0.5*position*position/t);
	g = g/(sqrt(2*MAT_PI*t));
	
	return g;	
}


float GaussKernelDerivative(int position , float t){
	float g;
	
	g = -position*exp(-0.5*position*position/t)/t;
	g = g/(sqrt(2*MAT_PI*t));
	
	return g;	
}


void ConvolveGauss(int type, float M, float t, float * in, float * out, int size){
	if(type==0){	
		for(int i=0;i<size;i++){
			if(i>=M && i<size-M){
			for(int k=-M;k<=M;k++){
					out[i]  += in[i-k]*GaussKernel(k,t);
				}
			}
		}					
	}
	if(type==1){
		for(int i=0;i<size;i++){
			if(i>=M && i<size-M){
			for(int k=-M;k<=M;k++){
					out[i] += in[i-k]*GaussKernelDerivative(k,t);
				}
			}
			out[i] = out[i] * size;
		}							
	}
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


void ExtremaFinding(float *vec_in, int N, float sound, float dt,float & sat, float & tau , float & error, std::string extremafile){		
	ofstream data_extrema;
	data_extrema.open(extremafile);
	
	data_extrema << "#pos_max" << "\t" << "Max" <<"\t"<< "pos_min" <<"\t"<< "Min"<<endl;	 
	
	int W;
	int pos_max, pos_min;
	float maximum,minimum;
	
	W = floor( 1.2*2*MAT_PI/(RealFreq(sound, 1.0, 1.0, 1)*dt));	
	int k = 0;
	int M = ceil(0.5*dt*N*RealFreq(sound, 1.0, 1.0, 1)/MAT_PI);
	float *vec_max;			
	vec_max =(float*) calloc (M,sizeof(float));
	float *vec_pos;			
	vec_pos =(float*) calloc (M,sizeof(float));
	
	for(int shift=0; shift < N-W ; shift += W){
		maximum =  *max_element( vec_in + shift ,vec_in + shift + W );
		pos_max =   max_element( vec_in + shift ,vec_in + shift + W ) - vec_in; 
		
		minimum =  *min_element( vec_in + shift ,vec_in + shift + W );
		pos_min =   min_element( vec_in + shift ,vec_in + shift + W ) - vec_in; 
		
		data_extrema << pos_max*dt << "\t" << maximum <<"\t"<< pos_min*dt <<"\t"<< minimum  <<endl;	 	
		vec_max[k] = maximum;
		vec_pos[k] = pos_max*dt;
		k++;
	}
	sat =  *max_element( vec_max ,vec_max + M );	
	for(int i=1;i<M-1;i++){
		if( vec_max[i] > 0.99*sat ){
			tau  = 	vec_pos[i];
			error = 0.5*(vec_pos[i+1]-vec_pos[i]);
			break;
		}
	}
	data_extrema.close();		
}

void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile){
	ofstream data_shock;
	data_shock.open(shockfile);
	for ( int i = 0; i < N; i++ )
	{
		if(i>=1){
			float schockD=0.0;
			schockD = (in[i]-in[i-1])/dx;
			if(abs(schockD) >=20){
				data_shock  <<t<<"\t"<< i*dx <<endl;	
			}
		}
	}
	data_shock.close();	
}

void Autocorrelation(float * out_gamma ,float * in , int crop, int size){
	int M = size - crop;
	float in_crop[M];
	for(int k =0;k < M;k++){
		in_crop[k] = in[k+crop];		
	}
	float sum;
	for(int lag=0; lag < M ;lag++){		
		for(int t=0; t < M ;t++){
			sum += in_crop[ t ] * in_crop[ (t + lag)%M];
		}
		out_gamma[lag] = sum;
		sum=0.0;
	}
}

void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out ){
	//second order method
	//f[tempo][posicao]
	for(int i=1;i<size_rows-1;i++){
		for(int j=0;j<size_cols;j++){
			df_out[i][j] = (-0.5*f_in[i-1][j]+0.5*f_in[i+1][j])/dt;
		}
	}

	for(int j=0;j<size_cols;j++){
		df_out[0][j]           = (-1.5*f_in[0][j]+2.0*f_in[1][j]-0.5*f_in[2][j])/dt;
		df_out[size_rows-1][j] = ( 0.5*f_in[size_rows-1-2][j]-2.0*f_in[size_rows-1-1][j]+1.5*f_in[size_rows-1][j])/dt;
	}
}


     
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out ){
	//second order method
	//f[tempo][posicao]
	for(int i=0;i<size_rows;i++){
		for(int j=1;j<size_cols-1;j++){
			df_out[i][j] = (-0.5*f_in[i][j-1]+0.5*f_in[i][j+1])/dt;
		}
	}
	for(int i=0;i<size_rows;i++){
		df_out[i][0]           = (-1.5*f_in[i][0]+2.0*f_in[i][1]-0.5*f_in[i][2])/dt;
		df_out[i][size_cols-1] = ( 0.5*f_in[i][size_cols-1-2]-2.0*f_in[i][size_cols-1-1]+1.5*f_in[i][size_cols-1])/dt;
	}
}



float CartDistance(float x , float y , float z, float X , float Y , float Z ){
	float rsquared;
	rsquared = pow(X-x,2) + pow(Y-y,2) + pow(Z-z,2);
	return sqrt(rsquared);
}

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
	float q =-1.0;
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
