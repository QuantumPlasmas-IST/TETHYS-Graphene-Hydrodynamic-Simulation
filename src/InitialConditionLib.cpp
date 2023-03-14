/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "InitialConditionLib.h"

void InitialCondition::Rand(Fluid1D &fluid_class) {

}

void InitialCondition::Rand(Fluid2D &fluid_class) {
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < fluid_class.Nx*fluid_class.Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ;
		fluid_class.Umain[c].n() = 1.0f + 0.005f * (noise - 0.5f);
		noise =  (float) rd()/maxrand ;
		fluid_class.Umain[c].tmp()=  .2f + 0.005f * (noise - 0.5f);
	}
}

void InitialCondition::Rand(DiracGraphene2D &fluid_class) {
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();
	for (int c = 0; c < fluid_class.Nx*fluid_class.Ny; c++ ){
		float noise;
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		//Den[c] = 1.0f + 0.05f * (noise - 0.5f);
		fluid_class.Umain[c].n()=1.0f + 0.05f * (noise - 0.5f);
		noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
		//HDen[c] = 1.0f + 0.05f * (noise - 0.5f);
		fluid_class.HoleUmain[c].n()=1.0f + 0.05f * (noise - 0.5f);
	}
}

void InitialCondition::Test(Fluid2D &fluid_class) {
	for (int i = 0; i < fluid_class.Nx; i++ ){
		for (int j=0; j<fluid_class.Ny; j++){
			float densi;
			if(i>=80&&i<=120&&j>=80&&j<=120){
				densi=0.2f;
			}
			else{
				densi=0.0f;
			}
			fluid_class.Umain[i + j * fluid_class.Nx].n() = 1.0f + densi;
			fluid_class.Umain[i + j * fluid_class.Nx].px() = 0.1f;
		}
	}
}

void InitialCondition::Wave(Fluid2D &fluid_class) {
	for (int i = 0; i < fluid_class.Nx; i++ ){
		for (int j=0; j<fluid_class.Ny; j++){
			fluid_class.Umain[i + j * fluid_class.Nx].n() = 1.0f+0.05f*sin(5.0f*2.0f*MAT_PI*i*fluid_class.dx/fluid_class.lengX) ;
			fluid_class.Umain[i + j * fluid_class.Nx].px() = 0.0f;
		}
	}
}

void InitialCondition::InitialCondGeneral(Fluid2D &fluid_class, function<float(float, float)> fden,function<float(float, float)> fvx, function<float(float, float)> fvy) {
	float x,y;
	for (int i = 0; i < fluid_class.Nx; i++ ){
		for (int j=0; j<fluid_class.Ny; j++){
			x=i*fluid_class.dx;
			y=j*fluid_class.dy;
			fluid_class.Umain[i + j * fluid_class.Nx].n() = fden(x,y);
			fluid_class.Umain[i + j * fluid_class.Nx].px() = fvx(x,y);
			fluid_class.Umain[i + j * fluid_class.Nx].py() = fvy(x,y);
		}
	}
}
