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


//<<<<<<< add_source_term
float DensityFlux(float den,float vel,float vel_snd,float vel_fer);

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer);

float DensitySource(float den,float vel,float vel_snd,float vel_fer);

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float mfp);



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
	float vel_snd;						    // Sound speed
	float vel_fer;							// Fermi velocity
	float mfp; 								// mean free path in units of GFET length
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
	
	if(argc==5){
		vel_snd = atof(argv[1]);
		vel_fer = atof(argv[2]);
		mfp     = atof(argv[3]); 
		data_save_mode = atoi(argv[4]);	// full data or light save option
		}
	else{
		cout << "Define S value: ";
		cin >> vel_snd;
		cout << "Define Vf value: ";
		cin >> vel_fer;
		cout << "Define mean free path value: ";
		cin >> mfp;
		cout << "Define data_save_mode value (0-> light save | 1-> full data):";
		cin >> data_save_mode;
		}
	
	
	// NEEDS REVIEWING
	/*......CFL routine to determine dt...............................*/	
	dx = leng / ( float ) ( Nx - 1 );
	if(vel_snd<5){
		dt = dx / (5*vel_snd);
	}
	else{
		if(vel_snd>8 && vel_snd<10){
			dt = dx / (30+3*vel_snd);
		}
		else{
			dt = dx / (20+2*vel_snd);		
			//dt = dx / (5+1.5*vel_snd);		
		}
	}
	/*................................................................*/
	
	
	/*.........Fixed or variable vel_snd value..............................*/
	float *arr_snd;							
	arr_snd =(float*) calloc (Nx,sizeof(float));	
	for(int i = 0; i<Nx  ;i++){
		//arr_snd[i]= vel_snd - 0.15*vel_snd*( dx*i- floor(dx*i) );
		arr_snd[i]=vel_snd;
	}
	/*................................................................*/


	/*.........Output files and streams...............................*/
	
	// density(x,t)
	string densityfile = "density_" + to_string(vel_snd)+ ".dat" ;
	densityfile.erase (densityfile.end()-9, densityfile.end()-5);
	ofstream data_density;
	data_density.open (densityfile);
	data_density << fixed ;
	data_density << setprecision(6);
	// velocity(x,t)	
	string velocityfile = "velocity_" + to_string(vel_snd)+ ".dat" ;
	velocityfile.erase (velocityfile.end()-9, velocityfile.end()-5);
	ofstream data_velocity;
	data_velocity.open (velocityfile);
	data_velocity << fixed ;
	data_velocity << setprecision(6);
	// current(x,t)	
	string currentfile = "current_" + to_string(vel_snd)+ ".dat" ;
	currentfile.erase (currentfile.end()-9, currentfile.end()-5);
	ofstream data_current;
	data_current.open (currentfile);
	data_current << fixed ;
	data_current << setprecision(6);	

	// time density(L,t)-1=U(L,t) current(0,t) electric_dipole_moment(t)  derivative_electric_dipole_moment(t)
	string electrofile = "electro_" + to_string(vel_snd)+ ".dat" ;
	electrofile.erase (electrofile.end()-9, electrofile.end()-5);
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific; 
	// time density(L,t) velocity(L,t) density(0,t) velocity(0,t)
	string slicefile = "slice_" + to_string(vel_snd)+ ".dat" ;
	slicefile.erase (slicefile.end()-9, slicefile.end()-5);							
	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	/*................................................................*/

	
	
//	cout << "\n*******************************************************"<< endl;
	cout << "Sound speed S\t"<< vel_snd <<endl;
	cout <<"dt= "<<dt<<"\tdx= "<<dx<<endl;
	cout << "Predicted w'= "<< RealFreq(vel_snd,1.0,1.0,1) << "\t1/w'= "<< 1.0/RealFreq(vel_snd,1.0,1.0,1)  << endl;
	cout << "Predicted w''= "<< ImagFreq(vel_snd,1.0,1.0) <<"\t1/w''= "<< 1.0/ImagFreq(vel_snd,1.0,1.0) <<endl;
	
	logfile << "#vel_snd \t dt \t dx \t w' \t w'' " << endl;
	logfile << vel_snd <<"\t"<< dt <<"\t"<< dx <<"\t"<< RealFreq(vel_snd,1.0,1.0,1) <<"\t"<< ImagFreq(vel_snd,1.0,1.0) ;
	
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
				- ( 0.5*dt/dx ) * ( DensityFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - DensityFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 
				+ ( 0.5*dt    ) * DensitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer) ;
			vel_mid[i] = 0.5*( vel[i] + vel[i+1] )
				- ( 0.5*dt/dx ) * ( VelocityFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - VelocityFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer, mfp) ;
		}
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - DensityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * DensitySource(den[i],vel[i],arr_snd[i], vel_fer);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - VelocityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * VelocitySource(den[i],vel[i],arr_snd[i], vel_fer, mfp);
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

float DensityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f1;
	
	f1 = den*vel;
	
	return f1;
}


float VelocityFlux(float den,float vel,float vel_snd,float vel_fer){
	float f2;
	f2 = 0.25*vel*vel + vel_fer*vel_fer*0.5*log(den) + 2*vel_snd*vel_snd*sqrt(den); 
	return f2;
}

float DensitySource(float den,float vel,float vel_snd,float vel_fer){
	float Q1=0.0;
return Q1;	
}

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float mfp){
	float Q2=0.0;
	if(mfp==0.0){
		Q2=0.0;
	}
	else{
		Q2=-1.0*(vel-1)/mfp;
		}
return Q2;
}

