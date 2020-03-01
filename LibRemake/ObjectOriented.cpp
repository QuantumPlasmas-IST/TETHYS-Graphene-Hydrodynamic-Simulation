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

//#include <H5Cpp.h>


//#include "TestLib.h"
#include "TethysLib.h"

using namespace std;


int main(int argc, char **argv){

float T_max=10.0;
int data_save_mode=0;
float dx;								// spatial discretisation
float dt;	
float t=0.0;
int time_step=0;

	GrapheneFluid1D Grafeno(201);	
	Grafeno.CFLCondition();
	Grafeno.SetVelFer(10.0);
	Grafeno.SetVelSnd(43.0);
	Grafeno.SetSound();
	
	dx=Grafeno.GetDx();
	dt=Grafeno.GetDt(); 
	cout << dx<<"\t"<< dt<<endl;
	
	/*.........Output files and streams...............................*/
	string str_snd = to_string(Grafeno.GetVelSnd());
	str_snd.erase(str_snd.end()-4,str_snd.end());
	string str_fer = to_string(Grafeno.GetVelFer());
	str_fer.erase(str_fer.end()-4,str_fer.end());
	string str_col_freq = to_string(Grafeno.GetColFreq());
	str_col_freq.erase(str_col_freq.end()-4,str_col_freq.end());
	string nam_post = "S="+str_snd+"vF="+str_fer+"l="+str_col_freq;
	// time density(L,t)-1=U(L,t) current(0,t) electric_dipole_moment(t)  derivative_electric_dipole_moment(t)
	//string electrofile = "electro_" + nam_post + ".dat" ;
	//ofstream data_electro;
	//data_electro.open (electrofile);
	//data_electro << scientific; 
	// time density(L,t) velocity(L,t) density(0,t) velocity(0,t)
	string slicefile = "slice_" + nam_post + ".dat" ;
	ofstream data_slice;
	data_slice.open (slicefile);
	data_slice << scientific; 
	/*................................................................*/


	
	Grafeno.InitialCondRand();
	
	
	
	


	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	
	
	while(t<=T_max && isfinite(Grafeno.vel[(Grafeno.SizeX()-1)/2]))
	{	
		++time_step;
		t += dt;
		
		Grafeno.Richtmyer();
		
		// Impose boundary conditions
		Grafeno.BoundaryCond(3);
		
		// Applying average filters for smoothing 	
		Grafeno.Smooth(2);
		
		/*
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
		}*/
		//Record end points
		data_slice <<t<<"\t"<< Grafeno.den_cor[Grafeno.SizeX()-1] <<"\t"<< Grafeno.vel_cor[Grafeno.SizeX()-1] <<"\t"<< Grafeno.den_cor[0] <<"\t" << Grafeno.vel_cor[0] <<"\n";
		//Record electric quantities
		//data_electro <<t<<"\t"<< graph.den_cor[Nx-1]-1.0 <<"\t"<< graph.den_cor[0]*graph.vel_cor[0] <<"\n";//<<"\t"<<  TotalElectricDipole(Nx,dx,den_cor)<<"\t"<<  DtElectricDipole(Nx,dx,cur_cor) <<"\t"<< KineticEnergy(Nx,dx, den_cor, vel_cor)  <<"\n";
	}
	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";

return 0;	
}


