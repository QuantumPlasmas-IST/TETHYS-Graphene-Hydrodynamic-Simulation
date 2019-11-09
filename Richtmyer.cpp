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

#include "dyakonovshur.h"


using namespace std;


float DensityFlux(float den,float vel,float s,float vF);

float VelocityFlux(float den,float vel,float s,float vF);

float DensitySource(float den,float vel,float s,float vF);

float VelocitySource(float den,float vel,float s,float vF);




int main(int argc, char **argv){

	cout << "*******************************************************"<<endl;
	cout << "************** ELECTRON FLOW SIMULATION ***************"<<endl;
	cout << "**************     Richtmyer method     ***************"<<endl;
	cout << "**************      Version 1.0.0       ***************"<<endl;
	cout << "*******************************************************"<<endl;


	/*......TIME stamp for the logfile................................*/
	ofstream logfile;
	logfile.open("Simulation.log",std::ios_base::app);
	time_t time_raw;
	struct tm * time_info;
	time (&time_raw);
	time_info = localtime (&time_raw);
	logfile << "\n#Simulation @ " << asctime(time_info) ;
	/*................................................................*/
	
	int Nx=201; 							// number of spatial points
	float t=0.0,leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
	float S;							    // Sound speed

	float *den;							 	//density field
	den =(float*) calloc (Nx,sizeof(float));
	float *den_mid;							//density auxiliary vector for midpoint calculation 
	den_mid = (float*) calloc (Nx-1,sizeof(float));
	float *vel;								//velocity field
 	vel = (float*) calloc (Nx,sizeof(float));
	float *vel_mid;							//velocity auxiliary vector for midpoint calculation 
	vel_mid = (float*) calloc (Nx-1,sizeof(float));
	
	float *den_cor;							//density corrected after average filter 
	den_cor = (float*) calloc (Nx,sizeof(float));
	float *vel_cor;							//velocity corrected after average filter 
	vel_cor = (float*) calloc (Nx,sizeof(float));
 	float *cur_cor;							//current density (n*v) corrected after average filter 
	cur_cor = (float*) calloc (Nx,sizeof(float));
 	
 	
 	int data_save_mode=0;
	
	if(argc!=1){
		
		S = atof(argv[1]);
		data_save_mode = atoi(argv[2]);	// full data or light save option
	
		}
	else{
		cout << "Define S value: ";
		cin >> S;
		}
	
	
	// NEEDS REVIEWING
	/*......CFL routine to determine dt...............................*/	
	dx = leng / ( float ) ( Nx - 1 );
	if(S<5){
		dt = dx / (5*S);
	}
	else{
		if(S>8 && S<10){
			dt = dx / (30+3*S);
		}
		else{
			dt = dx / (20+2*S);		
			//dt = dx / (5+1.5*S);		
		}
	}
	/*................................................................*/
	
	
	/*.........Fixed or variable S value..............................*/
	float *s;							
	s =(float*) calloc (Nx,sizeof(float));	
	for(int i = 0; i<Nx  ;i++){
		//s[i]= S - 0.15*S*( dx*i- floor(dx*i) );
		s[i]=S;
	}
	/*................................................................*/


	/*.........Output files and streams...............................*/
	
	// density(x,t)
	string densityfile = "density_" + to_string(S)+ ".dat" ;
	densityfile.erase (densityfile.end()-9, densityfile.end()-5);
	ofstream data_density;
	data_density.open (densityfile);
	data_density << fixed ;
	data_density << setprecision(6);
	// velocity(x,t)	
	string velocityfile = "velocity_" + to_string(S)+ ".dat" ;
	velocityfile.erase (velocityfile.end()-9, velocityfile.end()-5);
	ofstream data_velocity;
	data_velocity.open (velocityfile);
	data_velocity << fixed ;
	data_velocity << setprecision(6);
	// current(x,t)	
	string currentfile = "current_" + to_string(S)+ ".dat" ;
	currentfile.erase (currentfile.end()-9, currentfile.end()-5);
	ofstream data_current;
	data_current.open (currentfile);
	data_current << fixed ;
	data_current << setprecision(6);	

	// time density(L,t)-1=U(L,t) current(0,t) electric_dipole_moment(t)  derivative_electric_dipole_moment(t)
	string electrofile = "electro_" + to_string(S)+ ".dat" ;
	electrofile.erase (electrofile.end()-9, electrofile.end()-5);
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific; 
	// time density(L,t) velocity(L,t) density(0,t) velocity(0,t)
	string slicefile = "slice_" + to_string(S)+ ".dat" ;
	slicefile.erase (slicefile.end()-9, slicefile.end()-5);							
	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	/*................................................................*/

	
	
//	cout << "\n*******************************************************"<< endl;
	cout << "Sound speed S\t"<< S <<endl;
	cout <<"dt= "<<dt<<"\tdx= "<<dx<<endl;
	cout << "Predicted w'= "<< RealFreq(S,1.0,1.0,1) << "\t1/w'= "<< 1.0/RealFreq(S,1.0,1.0,1)  << endl;
	cout << "Predicted w''= "<< ImagFreq(S,1.0,1.0) <<"\t1/w''= "<< 1.0/ImagFreq(S,1.0,1.0) <<endl;
	
	logfile << "#S \t dt \t dx \t w' \t w'' " << endl;
	logfile << S <<"\t"<< dt <<"\t"<< dx <<"\t"<< RealFreq(S,1.0,1.0,1) <<"\t"<< ImagFreq(S,1.0,1.0) ;
	
//	cout << "*******************************************************"<< endl;
	
	float T_max=10.0;
	
	cout <<"Determined maximum simulated time\t" <<T_max<<endl;
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	InitialCondRand(Nx, dx, den, vel);
	BoundaryCond(3, Nx, den, vel);
	////////////////////////////////////////////////////////////////////
	
	if(data_save_mode){
		for(int i = 0; i<Nx  ;i++)
		{
			data_density   <<  den[i] <<"\t";
			data_current   <<  vel[i]*den[i] <<"\t";
			data_velocity  <<  vel[i] <<"\t";
		}
	}
	
	cout << "Running"<<endl;
	
	int time_step=0;
	
	while(t<=T_max && isfinite(vel[(Nx-1)/2]))
	{	
		++time_step;
		t += dt;
		
		if(data_save_mode && time_step % 35 == 0 ){
			data_density  << "\n";
			data_current  << "\n";
			data_velocity << "\n";
		}
		
		//
		//  Half step calculate density and velocity at time k+0.5 at the spatial midpoints
		//
		for ( int i = 0; i < Nx - 1; i++ )
		{
			den_mid[i] = 0.5*( den[i] + den[i+1] )
				- ( 0.5*dt/dx ) * ( DensityFlux(den[i+1],vel[i+1],s[i+1],10.0) - DensityFlux(den[i],vel[i],s[i],10.0) ) 
				+ ( 0.5*dt    ) * DensitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),s[i],10.0) ;
			vel_mid[i] = 0.5*( vel[i] + vel[i+1] )
				- ( 0.5*dt/dx ) * ( VelocityFlux(den[i+1],vel[i+1],s[i+1],10.0) - VelocityFlux(den[i],vel[i],s[i],10.0) ) 
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),s[i],10.0) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],s[i],10.0) - DensityFlux(den_mid[i-1],vel_mid[i-1],s[i-1],10.0) )
							+  dt * DensitySource(den[i],vel[i],s[i],10.0);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],s[i],10.0) - VelocityFlux(den_mid[i-1],vel_mid[i-1],s[i-1],10.0) )
							+  dt * VelocitySource(den[i],vel[i],s[i],10.0);
		}
		
		// Impose boundary conditions
		BoundaryCond(3, Nx, den, vel);
		
		// Applying average filters for smoothing 
		AverageFilter( den ,den_cor, Nx , 2);	
		AverageFilter( vel ,vel_cor, Nx , 2);
		
		// calculate current density already smoothed			
		for ( int i = 0; i < Nx; i++ )
		{	
			cur_cor[i] = 	vel_cor[i]*den_cor[i];	
			if(data_save_mode){
			//Record full data
				if(time_step % 35 == 0){
					data_density   << den_cor[i] <<"\t";
					data_current   << cur_cor[i] <<"\t";
					data_velocity  << vel_cor[i] <<"\t";
				}
			}
		}
	
		//Record end points
		data_slice <<t<<"\t"<< den_cor[Nx-1] <<"\t"<< vel_cor[Nx-1] <<"\t"<< den_cor[0] <<"\t" << vel_cor[0] <<"\n";
		//Record electric quantities
		data_electro <<t<<"\t"<< den_cor[Nx-1]-1.0 <<"\t"<< den_cor[0]*vel_cor[0] <<"\t"<<  TotalElectricDipole(Nx,dx,den_cor)<<"\t"<<  DtElectricDipole(Nx,dx,cur_cor) <<"\t"<< KineticEnergy(Nx,dx, den_cor, vel_cor)  <<"\n";
	}
	cout << "DONE!" <<endl;
	cout << "*******************************************************"<<endl;

	free(den);
	free(den_mid);
	free(den_cor);
	free(vel);
	free(vel_mid);
	free(vel_cor);
	free(cur_cor);

	data_density.close();
	data_velocity.close();
	data_current.close();
	data_slice.close();
	data_electro.close();
	
	
	
	
	return 0;
}

float DensityFlux(float den,float vel,float sound, float fermi){
	float f1;
	
	f1 = den*vel;
	
	return f1;
}

float VelocityFlux(float den,float vel,float sound, float fermi){
	float f2;
	
	f2 = 0.25*vel*vel + fermi*fermi*0.5*log(den) + 2*sound*sound*sqrt(den); 
	
	
	return f2;
}



float DensitySource(float den,float vel,float s,float vF){
	float Q1=0.0;
return Q1;	
}

float VelocitySource(float den,float vel,float s,float vF){
//	float Q2=0.0;
float tau = 2.0;
	float Q2=-1.0*(vel-1)/tau;
return Q2;
}

