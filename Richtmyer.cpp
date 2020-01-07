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

#include "H5Cpp.h"

#include "dyakonovshur.h"

using namespace H5;
using namespace std;
const H5std_string   FILE_NAME( "Richtmyer.h5" );
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);


float DensityFlux(float den,float vel,float vel_snd,float vel_fer);

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer);

float EnergyFlux(float den,float vel,float vel_snd,float vel_fer);

float DensitySource(float den,float vel,float vel_snd,float vel_fer);

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq);

float EnergySource(float den,float den_der,float vel,float vel_snd,float vel_fer);

int main(int argc, char **argv){
	/* Display name and version  */
    BannerDisplay();

	const int Nx=201; 							// number of spatial points
	float t=0.0;
	float T_max=10.0;
	const float leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
	float vel_snd;						    // Sound speed
	float vel_fer;							// Fermi velocity
	float col_freq; 								// mean free path in units of GFET length
	
	float *den;							 	//density field
	den =(float*) calloc (Nx,sizeof(float));
	float *den_mid;							//density auxiliary vector for midpoint calculation 
	den_mid = (float*) calloc (Nx-1,sizeof(float));
//	float *eng;							 	//energy density field
//	eng =(float*) calloc (Nx,sizeof(float));
//	float *eng_mid;							//energy density auxiliary vector for midpoint calculation 
//	eng_mid = (float*) calloc (Nx-1,sizeof(float));
	float *vel;								//velocity field
 	vel = (float*) calloc (Nx,sizeof(float));
	float *vel_mid;							//velocity auxiliary vector for midpoint calculation 
	vel_mid = (float*) calloc (Nx-1,sizeof(float));

	float *den_cor;							//density corrected after average filter 
	den_cor = (float*) calloc (Nx,sizeof(float));
	float *vel_cor;							//velocity corrected after average filter 
	vel_cor = (float*) calloc (Nx,sizeof(float));
 //	float *eng_cor;							//energy density corrected after average filter 
//	eng_cor = (float*) calloc (Nx,sizeof(float));
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

	
	
	/*.............  HDF5 file .......................................*/  	
	/*
	 * Create a file.
	 */
	H5File* hdf5file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
	/*
	 * Create the groups ("folders") in the file
	 */
	Group* grp_dat = new Group( hdf5file->createGroup( "/Data" ));
	Group* grp_den = new Group( hdf5file->createGroup( "/Data/Density" ));
	Group* grp_vel = new Group( hdf5file->createGroup( "/Data/Velocity" ));
	Group* grp_cur = new Group( hdf5file->createGroup( "/Data/Current" ));
	/*
	 * Create attributes 
	 */	
	hsize_t dim_atr[1] = { 1 };
	// Create the data space for the attribute.
	DataSpace atr_dataspace = DataSpace (1, dim_atr );
	// Create a group attribute. 
	Attribute atr_vel_snd  = grp_dat->createAttribute( "S parameter", hdf5_float, atr_dataspace);
	Attribute atr_vel_fer  = grp_dat->createAttribute( "Fermi velocity", hdf5_float, atr_dataspace);
	Attribute atr_col_freq = grp_dat->createAttribute( "Collision frequency", hdf5_float, atr_dataspace);
	Attribute atr_dx = grp_dat->createAttribute( "Space discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_dt = grp_dat->createAttribute( "Time discretisation step", hdf5_float, atr_dataspace);
	Attribute atr_total_time = grp_dat->createAttribute( "Total simulation time", hdf5_float, atr_dataspace);
	Attribute atr_num_time_steps = grp_dat->createAttribute( "Number of time steps", hdf5_int, atr_dataspace);
	Attribute atr_num_space_points = grp_dat->createAttribute( "Number of spatial points", hdf5_int, atr_dataspace);
	// Write the attribute data. 
	atr_col_freq.write(hdf5_float, &col_freq);
	atr_vel_fer.write( hdf5_float, &vel_fer);
	atr_vel_snd.write( hdf5_float, &vel_snd);
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_total_time.write( hdf5_float, &T_max);
	atr_num_space_points.write( hdf5_int, &Nx);
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
	
	const int RANK = 1; // 1D simulation
	hsize_t     dim_snap[1];              // dataset dimensions
	dim_snap[0] = Nx;
	DataSpace dataspace_den( RANK, dim_snap );
	DataSpace dataspace_vel( RANK, dim_snap );
	DataSpace dataspace_cur( RANK, dim_snap );

	/*................................................................*/
	
	
	WellcomeScreen(vel_snd, vel_fer, col_freq, dt, dx, T_max);
	RecordLogFile(vel_snd, vel_fer, col_freq, dt, dx, T_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	InitialCondRand(Nx, dx, den, vel);
	BoundaryCond(3, Nx, den, vel);
	
	//for(int i = 0; i<Nx  ;i++)
	//{
	//	eng[i]=1.0;
	//}
	////////////////////////////////////////////////////////////////////
	
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	int time_step=0;
	
	while(t<=T_max && isfinite(vel[(Nx-1)/2]))
	{	
		++time_step;
		t += dt;
		
		
		
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
			//eng_mid[i] = 0.5*( eng[i] + eng[i+1] )
			//	- ( 0.5*dt/dx ) * ( EnergyFlux(den[i+1],vel[i+1],arr_snd[i], vel_fer) - EnergyFlux(den[i],vel[i],arr_snd[i], vel_fer) ) 	
			//	+ ( 0.5*dt    ) * EnergySource(0.5*(den[i]+den[i+1]),0.5*(-1.0*den[i]+den[i+1])/dx,0.5*(vel[i]+vel[i+1]),arr_snd[i], vel_fer);
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
			//eng[i] = eng[i] - (dt/dx) * ( EnergyFlux(den_mid[i],vel_mid[i],arr_snd[i], vel_fer) - EnergyFlux(den_mid[i-1],vel_mid[i-1],arr_snd[i], vel_fer) )				
			//				+  dt * EnergySource(den[i],0.5*(-1.0*den_mid[i-1]+den_mid[i])/dx,vel[i],arr_snd[i], vel_fer);	
		}
		
		// Impose boundary conditions
		BoundaryCond(3, Nx, den, vel);
		
		// Applying average filters for smoothing 
		AverageFilter( den ,den_cor, Nx , 2);	
		AverageFilter( vel ,vel_cor, Nx , 2);
		//AverageFilter( eng ,eng_cor, Nx , 2);
		
		// calculate current density already smoothed			
		for ( int i = 0; i < Nx; i++ )
		{	
			cur_cor[i] = 	vel_cor[i]*den_cor[i];	
		}	
		if(data_save_mode && time_step % 35 == 0 ){
		//Record full data
			string str_time = to_string(time_step/35);
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = grp_den->createDataSet( name_dataset , hdf5_float, dataspace_den );
			dataset_den.write( den_cor, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_vel = grp_vel->createDataSet( name_dataset , hdf5_float, dataspace_vel );
			dataset_vel.write( vel_cor, hdf5_float );
			dataset_vel.close();	
			
			DataSet dataset_cur = grp_cur->createDataSet( name_dataset , hdf5_float, dataspace_cur );
			dataset_cur.write( cur_cor, hdf5_float );
			dataset_cur.close();
		}
		//Record end points
		data_slice <<t<<"\t"<< den_cor[Nx-1] <<"\t"<< vel_cor[Nx-1] <<"\t"<< den_cor[0] <<"\t" << vel_cor[0] <<"\n";
		//Record electric quantities
		data_electro <<t<<"\t"<< den_cor[Nx-1]-1.0 <<"\t"<< den_cor[0]*vel_cor[0] <<"\t"<<  TotalElectricDipole(Nx,dx,den_cor)<<"\t"<<  DtElectricDipole(Nx,dx,cur_cor) <<"\t"<< KineticEnergy(Nx,dx, den_cor, vel_cor)  <<"\n";
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	atr_num_time_steps.write(hdf5_int, &time_step);
	atr_num_time_steps.close();
	
	free(den);
	free(den_mid);
	free(den_cor);
	free(vel);
	free(vel_mid);
	free(vel_cor);
	free(cur_cor);
	data_slice.close();
	data_electro.close();

	grp_dat->close(); 
	grp_den->close(); 
	grp_vel->close();
	hdf5file->close();
	
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

