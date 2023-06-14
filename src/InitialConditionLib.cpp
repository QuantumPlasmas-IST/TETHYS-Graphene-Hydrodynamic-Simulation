/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "InitialConditionLib.h"

void InitialCondition::Rand(Fluid1D &fluid_class) {
	random_device rd;
	float maxrand;
	maxrand = (float) random_device::max();

	for (int i = 0; i < fluid_class.Nx; i++ ){
		float noise = (float) rd()/ maxrand ;
		fluid_class.Umain[i].n()= 1.0f + 0.0001f * (noise - 0.5f);
		fluid_class.Umain[i].v()= 0.0f;
	}
	fluid_class.SetSound();
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
		noise =  (float) rd()/maxrand ;
		fluid_class.Umain[c].n()=1.0f + 0.05f * (noise - 0.5f);
		noise =  (float) rd()/maxrand ;
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

void InitialCondition::InitialCondGeneral(Fluid2D &fluid_class, function<float(float, float)> fden,function<float(float, float)> fvx, function<float(float, float)> fvy, Geometry *Geom) {
	float x,y;
	for (int i = 0; i < fluid_class.Nx; i++ ){
		for (int j=0; j<fluid_class.Ny; j++){
			if( (Geom->dominio.dom[i + j * fluid_class.Nx] == 1) || (Geom->fronteira.edg[i + j * fluid_class.Nx] == 1) ){
				x=i*fluid_class.dx;
				y=j*fluid_class.dy;
				fluid_class.Umain[i + j * fluid_class.Nx].n() = fden(x,y);
				fluid_class.Umain[i + j * fluid_class.Nx].px() = fvx(x,y);
				fluid_class.Umain[i + j * fluid_class.Nx].py() = fvy(x,y);
			}else{
				x=i*fluid_class.dx;
				y=j*fluid_class.dy;
				fluid_class.Umain[i + j * fluid_class.Nx].n() = fden(x,y);
				fluid_class.Umain[i + j * fluid_class.Nx].px() = 0;
				fluid_class.Umain[i + j * fluid_class.Nx].py() = 0;
			}
		}
	}
}

void InitialCondition::InitialCondPointDen(DiracGraphene2D &fluid_class){
	for (int c = 0; c < fluid_class.Nx*fluid_class.Ny; c++ ){
		fluid_class.Umain[c].n()=1.0f;
		fluid_class.HoleUmain[c].n()=1.0f;
	}
	for (int c = 3*fluid_class.Ny/7; c < 4*fluid_class.Ny/7; c++){
		for (int g = 2*fluid_class.Nx/7; g < 3*fluid_class.Nx/7; g++){
			fluid_class.Umain[c*fluid_class.Nx+g].n()=1.2f;
			fluid_class.HoleUmain[c*fluid_class.Nx+g].n()=0.8f;
		}
	}
	for (int c = 3*fluid_class.Ny/7; c < 4*fluid_class.Ny/7; c++){
		for (int g = 4*fluid_class.Nx/7; g < 5*fluid_class.Nx/7; g++){
			fluid_class.Umain[c*fluid_class.Nx+g].n()=0.8f;
			fluid_class.HoleUmain[c*fluid_class.Nx+g].n()=1.2f;
		}
	}
}

void InitialCondition::InitialCondUniform(DiracGraphene2D &fluid_class){
	for (int c = 0; c < fluid_class.Nx*fluid_class.Ny; c++ ){
		fluid_class.Umain[c].n()=1.0f;
		fluid_class.HoleUmain[c].n()=1.0f;
	}
}

void InitialCondition::InitialCondUniform(DiracGraphene2D &fluid_class, Geometry *Geom){
	for (int c = 0; c < fluid_class.Nx*fluid_class.Ny; c++ ){
		if(Geom->dominio.dom[c] == 1 || Geom->fronteira.edg[c] == 1){
			fluid_class.Umain[c].n()=1.0f;
			fluid_class.HoleUmain[c].n()=1.0f;
		}
	}
}

void InitialCondition::InitialCondTest(Fluid1D &fluid_class) {
	for (int i = 0; i < fluid_class.Nx; i++ ){
		fluid_class.Umain[i].v()= 1.0f/(1.0f+5.0f* pow(cosh((i*fluid_class.dx-0.5f)*12.0f),2.f));
		fluid_class.Umain[i].n()= 0.2f+0.2f/ pow(cosh((i*fluid_class.dx-0.5f)*12.0f),2.f); //(i>3*Nx/8 && i<5*Nx/8 ) ? 1.0f : 0.1f; //0.2f+0.2f/ pow(cosh((i*dx-0.5f)*12.0f),2);//
	}
	fluid_class.SetSound();
}

void InitialCondition::InitialCondGeneral(Fluid1D &fluid_class, function<float(float)> fden, function<float(float)> fvx) {
	float x;
	for (int i = 0; i < fluid_class.Nx; ++i) {
		x=i*fluid_class.dx;
		fluid_class.Umain[i].n()=fden(x);
		fluid_class.Umain[i].v()=fvx(x);
	}
	fluid_class.SetSound();
}
