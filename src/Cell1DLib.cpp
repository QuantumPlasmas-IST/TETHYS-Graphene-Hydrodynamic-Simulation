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

float CellHandler1D::Roe(int i) {
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


StateVec CellHandler1D::WENO3(int Nx, char side, char edge) { // I think it works ( o﹏o)

	// reconstructed value
	StateVec Uweno{};

	// interpolation vectors
	StateVec Uaux0{};
	StateVec Uaux1{};

	// smoothness parameters
	StateVec Ubeta0{};
	StateVec Ubeta1{};

	// weights
	StateVec Uomega0{};
	StateVec Uomega1{};

	StateVec Ualpha0{};
	StateVec Ualpha1{};

	// small positive number
	float eps = 1.0E-6;

	float d0 = 1.0f / 3.0f;
	float d1 = 2.0f / 3.0f;
	
	// ternary conditions will be used from now on to solve out of scope problems (this solution only works with periodic boundary condition)
	Ubeta0 = (index!=Nx-1) ? (U_ptr[index+1] - U_ptr[index]) * (U_ptr[index+1] - U_ptr[index]) : (U_ptr[0] - U_ptr[Nx-1]) * (U_ptr[0] - U_ptr[Nx-1]);
	Ubeta1 = (index!=0)    ? (U_ptr[index] - U_ptr[index-1]) * (U_ptr[index] - U_ptr[index-1]) : (U_ptr[0] - U_ptr[Nx-1]) * (U_ptr[0] - U_ptr[Nx-1]);

	switch(side) {
		case 'W' : // i - 1/2
			Uaux0 = (index!=Nx-1) ? 0.5f * (3.0f*U_ptr[index] - U_ptr[index+1]) : 0.5f * (3.0f*U_ptr[Nx-1] - U_ptr[0]);
			Uaux1 = 0.5f * (U_ptr[index] + U_ptr[index-1]);

			Ualpha0.n() = d0 / ((eps + Ubeta0.n()) * (eps + Ubeta0.n()));
			Ualpha0.v() = d0 / ((eps + Ubeta0.v()) * (eps + Ubeta0.v()));
			Ualpha1.n() = d1 / ((eps + Ubeta1.n()) * (eps + Ubeta1.n()));
			Ualpha1.v() = d1 / ((eps + Ubeta1.v()) * (eps + Ubeta1.v()));

			Uomega0 = Ualpha0 / (Ualpha0 + Ualpha1);
			Uomega1 = Ualpha1 / (Ualpha0 + Ualpha1);

			switch(edge) {
				case 'L' : { // i - 1
					// not the best solution performance-wise but it'll do for now
					CellHandler1D cellAux(index-1, fluid_ptr, U_ptr);
					Uweno = cellAux.WENO3(Nx,'E','L');
				}
				break;

				case 'R' : // i
					Uweno = Uomega0*Uaux0 + Uomega1*Uaux1;
				break;

				default: 0;
			}
		break;

		case 'E' : // i + 1/2
			Uaux0 = 0.5f * (U_ptr[index+1] + U_ptr[index]);
			Uaux1 = (index!=0) ? 0.5f * (3.0f*U_ptr[index] - U_ptr[index-1]) : 0.5f * (3.0f*U_ptr[0] - U_ptr[Nx-1]);

			Ualpha0.n() = d1 / ((eps + Ubeta0.n()) * (eps + Ubeta0.n()));
			Ualpha0.v() = d1 / ((eps + Ubeta0.v()) * (eps + Ubeta0.v()));
			Ualpha1.n() = d0 / ((eps + Ubeta1.n()) * (eps + Ubeta1.n()));
			Ualpha1.v() = d0 / ((eps + Ubeta1.v()) * (eps + Ubeta1.v()));

			Uomega0 = Ualpha0 / (Ualpha0 + Ualpha1);
			Uomega1 = Ualpha1 / (Ualpha0 + Ualpha1);

			switch(edge) {
				case 'L' : // i
					Uweno = Uomega0*Uaux0 + Uomega1*Uaux1;
				break;

				case 'R' : { // i + 1
					// not the best solution performance-wise but it'll do for now
					CellHandler1D cellAux(index+1, fluid_ptr, U_ptr);
					Uweno = cellAux.WENO3(Nx,'W','R');
				}
				break;

				default: 0;
			}
		break;

		default: 0;
	}

	// does not account for sound anisotropy
	Uweno.S() = fluid_ptr->GetVelSnd();
	return Uweno;
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

	Ureturn.n() = fluido->DensityFlux(R)  + fluido->DensityFlux(L);
	Ureturn.v() = fluido->VelocityFlux(R) + fluido->VelocityFlux(L);

	float A = max(fluido->JacobianSpectralRadius(L),fluido->JacobianSpectralRadius(R));
	Ureturn = Ureturn - A*(R-L);

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

StateVec NumericalSource::Average(Fluid1D *fluido, StateVec L, StateVec R) {
	StateVec Ureturn{};
	Ureturn.n()=fluido->DensitySource(0.5f*(L+R));
	Ureturn.v()=fluido->VelocitySource(0.5f*(L+R));
	return Ureturn;

}
