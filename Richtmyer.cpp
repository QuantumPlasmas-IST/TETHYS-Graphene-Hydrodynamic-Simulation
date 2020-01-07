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


float DensityFlux(float den,float vel,float vel_snd,float vel_fer);

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer);

float EnergyFlux(float den,float vel,float vel_snd,float vel_fer);

float DensitySource(float den,float vel,float vel_snd,float vel_fer);

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq);

float EnergySource(float den,float den_der,float vel,float vel_snd,float vel_fer);

int main(int argc, char **argv){
	/* Display name and version  */
    BannerDisplay();

	int Nx=201; 							// number of spatial points
	float t=0.0,leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
	float vel_snd;						    // Sound speed
	float vel_fer;							// Fermi velocity
	float col_freq; 								// mean free path in units of GFET length
	float *den;							 	//density field
	den =(float*) calloc (Nx,sizeof(float));
	float *den_mid;							//density auxiliary vector for midpoint calculation 
	den_mid = (float*) calloc (Nx-1,sizeof(float));
	float *eng;							 	//energy density field
	eng =(float*) calloc (Nx,sizeof(float));
	float *eng_mid;							//energy density auxiliary vector for midpoint calculation 
	eng_mid = (float*) calloc (Nx-1,sizeof(float));
	float *vel;								//velocity field
 	vel = (float*) calloc (Nx,sizeof(float));
	float *vel_mid;							//velocity auxiliary vector for midpoint calculation 
	vel_mid = (float*) calloc (Nx-1,sizeof(float));

	float *den_cor;							//density corrected after average filter 
	den_cor = (float*) calloc (Nx,sizeof(float));
	float *vel_cor;							//velocity corrected after average filter 
	vel_cor = (float*) calloc (Nx,sizeof(float));
 	float *eng_cor;							//energy density corrected after average filter 
	eng_cor = (float*) calloc (Nx,sizeof(float));
 	float *cur_cor;							//current density (n*v) corrected after average filter 
	cur_cor = (float*) calloc (Nx,sizeof(float));
 	
 	
 	
 	int data_save_mode=0;
	
	if(argc==5){
		vel_snd = atof(argv[1]);
		vel_fer = atof(argv[2]);
		col_freq = atof(argv[3]); 
		data_save_mode = atoi(argv[4]);	// full data or light save option
		}
	else{
		cout << "Define S value: ";
		cin >> vel_snd;
		cout << "Define vF value: ";
		cin >> vel_fer;
		cout << "Define collision frequency: ";
		cin >> col_freq;
		cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
		cin >> data_save_mode;
		}
	
	/*......CFL routine to determine dt...............................*/	
	dx = leng / ( float ) ( Nx - 1 );
	dt = TimeStepCFL(dx, vel_snd, vel_fer);
	
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
	
	string str_snd = to_string(vel_snd);
	str_snd.erase(str_snd.end()-4,str_snd.end());
	string str_fer = to_string(vel_fer);
	str_fer.erase(str_fer.end()-4,str_fer.end());
	string str_col_freq = to_string(col_freq);
	str_col_freq.erase(str_col_freq.end()-4,str_col_freq.end());
	
	string nam_post = "S="+str_snd+"vF="+str_fer+"l="+str_col_freq;
		
	// density(x,t)
	string densityfile = "density_" + nam_post + ".dat" ;
	ofstream data_density;
	// velocity(x,t)	
	string velocityfile = "velocity_" + nam_post + ".dat" ;
	ofstream data_velocity;
	// current(x,t)	
	string currentfile = "current_" + nam_post + ".dat" ;
	ofstream data_current;
	if(data_save_mode){
		data_density.open (densityfile);
		data_density << fixed ;
		data_density << setprecision(6);
		data_velocity.open (velocityfile);
		data_velocity << fixed ;
		data_velocity << setprecision(6);			
		data_current.open (currentfile);
		data_current << fixed ;
		data_current << setprecision(6);		
	}


	// time density(L,t)-1=U(L,t) current(0,t) electric_dipole_moment(t)  derivative_electric_dipole_moment(t)
	string electrofile = "electro_" + nam_post + ".dat" ;
	ofstream data_electro;
	data_electro.open (electrofile);
	data_electro << scientific; 
	// time density(L,t) velocity(L,t) density(0,t) velocity(0,t)
	string slicefile = "slice_" + nam_post + ".dat" ;
	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	/*................................................................*/

	
	float T_max=10.0;
	
	WellcomeScreen(vel_snd, vel_fer, col_freq, dt, dx, T_max);
	RecordLogFile(vel_snd, vel_fer, col_freq, dt, dx, T_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	InitialCondRand(Nx, dx, den, vel);
	BoundaryCond(3, Nx, den, vel);
	
	for(int i = 0; i<Nx  ;i++)
	{
		eng[i]=1.0;
	}
	////////////////////////////////////////////////////////////////////
	
	if(data_save_mode){
		for(int i = 0; i<Nx  ;i++)
		{
			data_density   <<  den[i] <<"\t";
			data_current   <<  vel[i]*den[i] <<"\t";
			data_velocity  <<  vel[i] <<"\t";
		}
	}
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
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
				+ ( 0.5*dt    ) * VelocitySource(0.5*(den[i]+den[i+1]),0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer, col_freq) ;
			/* NEW ENERGY FLUX */				
			eng_mid[i] = 0.5*( eng[i] + eng[i+1] )
				- ( 0.5*dt/dx ) * ( EnergyFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - EnergyFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 	
				+ ( 0.5*dt    ) * EnergySource(0.5*(den[i]+den[i+1]),0.5*(-1.0*den[i]+den[i+1])/dx,0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer);
		}
		
		
		//
		// Remaining step 
		//
		for ( int i = 1; i < Nx - 1; i++ )
		{
			den[i] = den[i] - (dt/dx) * ( DensityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - DensityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * DensitySource(den[i],vel[i],arr_snd[i], vel_fer);
			vel[i] = vel[i] - (dt/dx) * ( VelocityFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - VelocityFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )
							+  dt * VelocitySource(den[i],vel[i],arr_snd[i], vel_fer, col_freq);
			/* NEW ENERGY FLUX */				
			eng[i] = eng[i] - (dt/dx) * ( EnergyFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - EnergyFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )				
							+  dt * EnergySource(den[i],0.5*(-1.0*den_mid[i-1]+den_mid[i])/dx,vel[i],arr_snd[i], vel_fer);	
		}
		
		// Impose boundary conditions
		BoundaryCond(3, Nx, den, vel);
		
		// Applying average filters for smoothing 
		AverageFilter( den ,den_cor, Nx , 2);	
		AverageFilter( vel ,vel_cor, Nx , 2);
		AverageFilter( eng ,eng_cor, Nx , 2);
		
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
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	free(den);
	free(den_mid);
	free(den_cor);
	free(vel);
	free(vel_mid);
	free(vel_cor);
	free(cur_cor);

	if(data_save_mode){
		data_density.close();
		data_velocity.close();
		data_current.close();
	}	
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

float EnergyFlux(float den,float vel,float vel_snd,float vel_fer){
	float f3;
	f3 = vel*pow(den,1.5);
 	return f3;
}


float DensitySource(float den,float vel,float vel_snd,float vel_fer){
	float Q1=0.0;
return Q1;	
}

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq){
	float Q2=0.0;
	Q2=-1.0*col_freq*(vel-1);
return Q2;
}

float EnergySource(float den,float den_der,float vel,float vel_snd,float vel_fer){
	float Q3=0.0;
	Q3=pow(vel_snd/vel_fer,2)*den*vel*den_der;
return Q3;
}

