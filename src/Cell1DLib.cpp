//
// Created by pcosme on 04/02/2022.
//

#include "Cell1DLib.h"

CellHandler1D::CellHandler1D(int i, Fluid1D * ptr_fluid, StateVec * ptr_state) {
	index=i;
	fluid_ptr=ptr_fluid;
	U_ptr=ptr_state;
}

CellHandler1D::CellHandler1D(int i, StateVec * ptr_state) {
	index=i;
	fluid_ptr= nullptr;
	U_ptr=ptr_state;
}

float CellHandler1D::VanLeer(StateVec*Uin,int i) {  //TODO cada variável tem o seu limitador ATENÇAO
	float denom,numer,r,f;
	//float limit=2.0f;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];

	numer=Numer.v();
	denom=Denom.v();

	r=numer/(denom+tolerance);
	f=(r+abs(r))/(1+abs(r));
	return f;
}


float CellHandler1D::VanLeer(int i) {  //TODO cada variável tem o seu limitador ATENÇAO
	float denom,numer,r,f;
	//float limit=2.0f;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];

	numer=Numer.v();
	denom=Denom.v();

	r=numer/(denom+tolerance);
	f=(r+abs(r))/(1+abs(r));
	return f;
}

StateVec CellHandler1D::VanLeerU(int i) {
	StateVec Ureturn{};
	float denom,numer,r,f_vel=0,f_den=0;
	float limit=2.0f;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];

	/////////////////////////////VELOCITY
	numer=Numer.v();
	denom=Denom.v();
	if(denom!=0){
		r=numer/denom;
		f_vel=(r+abs(r))/(1+abs(r));
	} else{
		if(numer<=tolerance){
			f_vel=0;
		} else{
			f_vel=limit;
		}
	}
	/////////////////////////////Density
	numer=Numer.n();
	denom=Denom.n();
	if(denom!=0){
		r=numer/denom;
		f_den=(r+abs(r))/(1+abs(r));
	} else{
		if(numer<=tolerance){
			f_den=0;
		} else{
			f_den=limit;
		}
	}
	Ureturn.v()=f_vel;
	Ureturn.n()=f_den;
	return Ureturn;
}

float CellHandler1D::Roe(int i) { //TODO cut off
	float denom,numer,r,f;
	float limit=1.0f;
	float tolerance=1E-6;
	StateVec Numer = U_ptr[i]-U_ptr[i-1];
	StateVec Denom = U_ptr[i+1]-U_ptr[i];
	numer=Numer.v();
	denom=Denom.v();


	r=numer/(denom+tolerance);
	f=max(0.0f,min(1.0f,r));
	return f;
}

StateVec CellHandler1D::TVD(char side, char edge) {
	StateVec Utvd{};//(U_ptr[0]);
	switch(side) {
		case 'W' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[index]   + 0.5f * VanLeer(index) * (U_ptr[index + 1] - U_ptr[index]);
					break;
				case 'R' :
					Utvd=U_ptr[index+1]   - 0.5f * VanLeer(index+1) * (U_ptr[index + 2] - U_ptr[index+1]);
					break;
				default: 0;
			}
			break;
		case 'E' :
			switch(edge) {
				case 'L' :
					Utvd=U_ptr[index-1]   + 0.5f * VanLeer(index-1) * (U_ptr[index] - U_ptr[index-1]);
					break;
				case 'R' :
					Utvd=U_ptr[index]   - 0.5f * VanLeer(index) * (U_ptr[index + 1] - U_ptr[index]);
					break;
				default: 0;
			}
			break;
		default: 0;
	}
	return Utvd;
}


StateVec CellHandler1D::TVD(StateVec * Uin,int pos,char side, char edge) {
	StateVec Utvd{};//(Uin[0]);
	switch(side) {
		case 'W' :
			switch(edge) {
				case 'L' :
					Utvd=Uin[pos]   + 0.5f * VanLeer(Uin,pos) * (Uin[pos + 1] - Uin[pos]);
					break;
				case 'R' :
					Utvd=Uin[pos+1]   - 0.5f * VanLeer(Uin,pos+1) * (Uin[pos + 2] - Uin[pos+1]);
					break;
				default: 0;
			}
			break;
		case 'E' :
			switch(edge) {
				case 'L' :
					Utvd=Uin[pos-1]   + 0.5f * VanLeer(Uin,pos-1) * (Uin[pos] - Uin[pos-1]);
					break;
				case 'R' :
					Utvd=Uin[pos]   - 0.5f * VanLeer(Uin,pos) * (Uin[pos + 1] - Uin[pos]);
					break;
				default: 0;
			}
			break;
		default: 0;
	}
	return Utvd;
}

StateVec CellHandler1D::UNO(char side, char edge) {
	StateVec Uuno(U_ptr[0]);
	return Uuno;
}

StateVec NumericalFlux::Average(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn{};
	Ureturn.n()=fluido->DensityFlux(0.5f*(L+R));
	Ureturn.v()=fluido->VelocityFlux(0.5f*(L+R));
	//Ureturn.S()=fluido->GetVelSnd();
	return Ureturn;
}

StateVec NumericalFlux::Central(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn{};
	Ureturn.n()=fluido->DensityFlux(R)+fluido->DensityFlux(L);
	Ureturn.v()=fluido->VelocityFlux(R)+fluido->VelocityFlux(L);
	//Ureturn.S()=fluido->GetVelSnd();
	float A =max(fluido->JacobianSpectralRadius(R),fluido->JacobianSpectralRadius(L));
	Ureturn.n() = Ureturn.n() - A*(R-L).n();
	Ureturn.v() = Ureturn.v() - A*(R-L).v();
	return 0.5f*Ureturn;
}

StateVec NumericalFlux::Characteristic(Fluid1D *fluido,  StateVec L, StateVec R) {
	StateVec Ureturn(fluido->DensityFlux(R)+fluido->DensityFlux(L),fluido->VelocityFlux(R)+fluido->VelocityFlux(L));

	StateVec u{};
	u= 0.5f*(L+R);

	float Fden,Fvel;

	Fden = fluido->DensityFlux(R)-fluido->DensityFlux(L);
	Fvel = fluido->VelocityFlux(R)-fluido->VelocityFlux(L);

	Ureturn.n() = Ureturn.n() - ( fluido->JacobianSignum(u,"11")*Fden + fluido->JacobianSignum(u,"12")*Fvel );
	Ureturn.v() = Ureturn.v() - ( fluido->JacobianSignum(u,"21")*Fden + fluido->JacobianSignum(u,"22")*Fvel );

	return 0.5f*Ureturn;
}

