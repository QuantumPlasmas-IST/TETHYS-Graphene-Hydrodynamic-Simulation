#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>   
#include <cassert>

#include <H5Cpp.h>

#include "Tethys1DLib.h"


using namespace H5;
using namespace std;
const H5std_string   FILE_NAME( "Richtmyer_DATA_SET.h5" );
const FloatType      hdf5_float(PredType::NATIVE_FLOAT);
const IntType        hdf5_int(PredType::NATIVE_INT);




int main(int argc, char **argv){
	/* Display name and version  */
    BannerDisplay();


	const int Npoints=201; 							// number of spatial points
	float t=0.0;
	float T_max=10.0;
//	const float leng=1.0;					// time variable and spatial Length
	float dx;								// spatial discretisation
	float dt;								// time step
//	float vel_snd;						    // Sound speed
//	float vel_fer;							// Fermi velocity
//	float col_freq; 								// mean free path in units of GFET length



GrapheneFluid1D	graph(Npoints);

 	
 	int data_save_mode=0;
	float input_vel_snd,input_vel_fer,input_col_freq,input_kin_vis;
	if(argc==6){
		input_vel_snd = atof(argv[1]);
		input_vel_fer = atof(argv[2]);
		input_col_freq = atof(argv[3]);
		input_kin_vis = atof(argv[4]);
		assert(atoi(argv[5])==0 || atoi(argv[5])==1);
		data_save_mode = atoi(argv[5]);	// full data or light save option
		graph.SetVelSnd(input_vel_snd);
		graph.SetVelFer(input_vel_fer);
		graph.SetKinVis(input_kin_vis);
		graph.SetColFreq(input_col_freq);
		}
	else{
		cout << "Define S value: "; // throw exceptions if the velocities or frequency are negative or if S<Vf
		cin >> input_vel_snd;
		graph.SetVelSnd(input_vel_snd);
		cout << "Define vF value: ";
		cin >> input_vel_fer; 
		graph.SetVelFer(input_vel_fer);
		cout << "Define kinetic viscosity: ";
		cin >> input_kin_vis;
		graph.SetColFreq(input_kin_vis);
		cout << "Define collision frequency: ";
		cin >> input_col_freq;
		graph.SetColFreq(input_col_freq);
		cout << "Define data_save_mode value (0-> light save | 1-> full data): ";
		cin >> data_save_mode;
		}
	
	/*......CFL routine to determine dt...............................*/	
	graph.CFLCondition();
	dx=graph.GetDx();
	dt=graph.GetDt();
	graph.SetTmax(10.0);
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	graph.SetSound();
	/*................................................................*/


	/*.........Output files and streams...............................*/
	
	string str_snd = to_string(graph.GetVelSnd());
	str_snd.erase(str_snd.end()-4,str_snd.end());
	string str_fer = to_string(graph.GetVelFer());
	str_fer.erase(str_fer.end()-4,str_fer.end());
	string str_kin_vis = to_string(graph.GetKinVis());
	str_kin_vis.erase(str_kin_vis.end()-4,str_kin_vis.end());
	string str_col_freq = to_string(graph.GetColFreq());
	str_col_freq.erase(str_col_freq.end()-4,str_col_freq.end());
	
	string nam_post = "S="+str_snd+"vF="+str_fer+"vis="+str_kin_vis+"l="+str_col_freq;
		


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
	// throw HDF5 exceptions ???
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
	atr_col_freq.write(hdf5_float, &input_col_freq);
	atr_vel_fer.write( hdf5_float, &input_vel_fer);
	atr_vel_snd.write( hdf5_float, &input_vel_snd);
	atr_dx.write(hdf5_float, &dx);
	atr_dt.write( hdf5_float, &dt);
	atr_total_time.write( hdf5_float, &T_max);
	atr_num_space_points.write( hdf5_int, &Npoints);
	atr_col_freq.close();
	atr_vel_fer.close();
	atr_vel_snd.close();
	atr_dx.close();
	atr_dt.close();
	atr_total_time.close();
	atr_num_space_points.close();
	
	const int RANK = 1; // 1D simulation
	hsize_t     dim_snap[1];              // dataset dimensions
	dim_snap[0] = Npoints; 
	DataSpace dataspace_den( RANK, dim_snap );
	DataSpace dataspace_vel( RANK, dim_snap );
	DataSpace dataspace_cur( RANK, dim_snap );

	/*................................................................*/
	
	//WellcomeScreen(graph);
	WellcomeScreen(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(),graph.GetKinVis(), dt, dx, T_max);
	RecordLogFile(graph.GetVelSnd(), graph.GetVelFer(), graph.GetColFreq(), dt, dx, T_max);
	
	////////////////////////////////////////////////////////////////////
	// Initialization	
	graph.InitialCondRand();
	graph.BoundaryCond(3);
	////////////////////////////////////////////////////////////////////
	
	
	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	int time_step=0;
	
	while(t<=T_max && isfinite(graph.vel[(graph.SizeX()-1)/2])) // throw exception para nan / inf 
	{	
		++time_step;
		t += dt;
		
		graph.Richtmyer();
		
		// Impose boundary conditions
		graph.BoundaryCond(3);
		
		// Applying average filters for smoothing 	
		graph.Smooth(2);
		
		
		if(data_save_mode && time_step % 35 == 0 ){
		//Record full data
			string str_time = to_string(time_step/35);
			string name_dataset = "snapshot_"+str_time;
			
			DataSet dataset_den = grp_den->createDataSet( name_dataset , hdf5_float, dataspace_den );
			dataset_den.write( graph.den_cor, hdf5_float );
			dataset_den.close();
			
			DataSet dataset_vel = grp_vel->createDataSet( name_dataset , hdf5_float, dataspace_vel );
			dataset_vel.write( graph.vel_cor, hdf5_float );
			dataset_vel.close();	
			
			DataSet dataset_cur = grp_cur->createDataSet( name_dataset , hdf5_float, dataspace_cur );
			dataset_cur.write( graph.cur_cor, hdf5_float );
			dataset_cur.close();
		}
		//Record end points
		data_slice <<t<<"\t"<< graph.den_cor[Npoints-1] <<"\t"<< graph.vel_cor[Npoints-1] <<"\t"<< graph.den_cor[0] <<"\t" << graph.vel_cor[0] <<"\n";
		//Record electric quantities
		float Q_net, I_avg, P_ohm;
		Q_net = graph.NetCharge();
		I_avg = graph.AverageCurrent(); 
		P_ohm = graph.OhmPower();
		data_electro <<t<<"\t"<< Q_net<<"\t"<<I_avg<<"\t"<<Q_net*Q_net*0.5 <<"\t"<<P_ohm<<"\t"<<graph.ElectricDipole()<<"\t"<< graph.ElectricDipoleVariation() <<"\n";
		//data_electro <<t<<"\t"<< graph.den_cor[Nx-1]-1.0 <<"\t"<< graph.den_cor[0]*graph.vel_cor[0] <<"\n";//<<"\t"<<  TotalElectricDipole(Nx,dx,den_cor)<<"\t"<<  DtElectricDipole(Nx,dx,cur_cor) <<"\t"<< KineticEnergy(Nx,dx, den_cor, vel_cor)  <<"\n";
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout<<"═══════════════════════════════════════════════════════════════════════════" <<endl;

	atr_num_time_steps.write(hdf5_int, &time_step);
	atr_num_time_steps.close();

	data_slice.close();
	data_electro.close();

	grp_dat->close(); 
	grp_den->close(); 
	grp_vel->close();
	hdf5file->close();
	
	return 0;
}




