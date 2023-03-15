//
// Created by pcosme on 14/03/23.
//

#include "DiffOperatorLib.h"

void DiffOperator::VelocityGradient(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	VelocityXGradient_bulk(fluid_class,Uarray,size_x,size_y); //já paralelo
	VelocityYGradient_bulk(fluid_class,Uarray,size_x,size_y); //já paralelo

	VelocityXGradient_top(fluid_class,Uarray,size_x,size_y);
	VelocityXGradient_bottom(fluid_class,Uarray,size_x,size_y);
	VelocityXGradient_left(fluid_class,Uarray,size_x,size_y);
	VelocityXGradient_right(fluid_class,Uarray,size_x,size_y);
	VelocityXGradient_corners(fluid_class,Uarray,size_x,size_y);

	VelocityYGradient_top(fluid_class,Uarray,size_x,size_y);
	VelocityYGradient_bottom(fluid_class,Uarray,size_x,size_y);
	VelocityYGradient_left(fluid_class,Uarray,size_x,size_y);
	VelocityYGradient_right(fluid_class,Uarray,size_x,size_y);
	VelocityYGradient_corners(fluid_class,Uarray,size_x,size_y);
}



void DiffOperator::VelocityXGradient_bulk(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	float dx=fluid_class.dx,dy=fluid_class.dy;
#pragma omp parallel for default(none) shared(fluid_class,Uarray,stride,size_x,size_y,dx,dy)
	for (int kp = 1 + size_x; kp <= size_x * size_y - size_x - 2; kp++) {
		int N,S,E,W;
		if (kp % stride != stride - 1 && kp % stride != 0) {
			N=kp+stride;
			S=kp-stride;
			E=kp+1;
			W=kp-1;
			float mE = fluid_class.DensityToMass(Uarray[E].n());
			float mW = fluid_class.DensityToMass(Uarray[W].n());
			float mN = fluid_class.DensityToMass(Uarray[N].n());
			float mS = fluid_class.DensityToMass(Uarray[S].n());

			float vxE = Uarray[E].px()/mE;
			float vxW = Uarray[W].px()/mW;
			float vxN = Uarray[N].px()/mN;
			float vxS = Uarray[S].px()/mS;

			Uarray[kp].dxvx() = (vxE - vxW) / (2.0f * dx);
			Uarray[kp].dyvx() = (vxN - vxS) / (2.0f * dy);
		}
	}
}






void DiffOperator::VelocityXGradient_top(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int i = 1; i <= size_x - 2; i++) { // topo rede principal, ou seja j=(size_y - 1)
		int top = i + (size_y - 1) * stride;
		int southsouth = i + (size_y - 3) * stride;
		S=top-stride;
		E=top+1;
		W=top-1;
		float mC = fluid_class.DensityToMass(Uarray[top].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());
		float mSS = fluid_class.DensityToMass(Uarray[southsouth].n());

		float vxC = Uarray[top].px()/mC;
		float vxE = Uarray[E].px()/mE;
		float vxW = Uarray[W].px()/mW;
		float vxSS = Uarray[southsouth].px()/mSS;
		float vxS = Uarray[S].px()/mS;

		Uarray[top].dxvx() = (vxE - vxW) / (2.0f * dx);
		Uarray[top].dyvx() = (3.0f * vxC - 4.0f * vxS + vxSS) /(2.0f * dy); //backward finite difference
	}
}

void DiffOperator::VelocityXGradient_bottom(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int i = 1; i <= size_x - 2; i++) { // fundo rede principal, ou seja j=0
		int bottom = i; //i+0*nx
		int northnorth = i + 2 * stride;
		N=bottom+stride;
		E=bottom+1;
		W=bottom-1;
		float mC = fluid_class.DensityToMass(Uarray[bottom].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mNN = fluid_class.DensityToMass(Uarray[northnorth].n());

		float vxC = Uarray[bottom].px()/mC;
		float vxE = Uarray[E].px()/mE;
		float vxW = Uarray[W].px()/mW;
		float vxNN = Uarray[northnorth].px()/mNN;
		float vxN = Uarray[N].px()/mN;

		Uarray[bottom].dxvx() = (vxE - vxW) / (2.0f * dx);
		Uarray[bottom].dyvx() = (-3.0f * vxC + 4.0f * vxN - vxNN) /(2.0f * dy); //backward finite difference
	}
}

void DiffOperator::VelocityXGradient_left(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int j = 1; j <= size_y - 2; j++) { //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j * stride;
		int easteast = left + 2;
		N=left+stride;
		S=left-stride;
		E=left+1;

		float mC = fluid_class.DensityToMass(Uarray[left].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mEE = fluid_class.DensityToMass(Uarray[easteast].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());

		float vxC = Uarray[left].px()/mC;
		float vxE = Uarray[E].px()/mE;
		float vxEE = Uarray[easteast].px()/mEE;
		float vxN = Uarray[N].px()/mN;
		float vxS = Uarray[S].px()/mS;

		Uarray[left].dxvx() = (-3.0f * vxC+ 4.0f * vxE - vxEE) /(2.0f * dx); //forward difference
		Uarray[left].dyvx() = (vxN - vxS) / (2.0f * dy); //OK
	}
}

void DiffOperator::VelocityXGradient_right(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int j = 1; j <= size_y - 2; j++) { //lado direito da rede principal ou seja i=(size_x-1)
		int right = (size_x - 1) + j * stride;
		int westwest = right - 2;
		N=right+stride;
		S=right-stride;
		W=right-1;
		float mC = fluid_class.DensityToMass(Uarray[right].n());
		float mWW = fluid_class.DensityToMass(Uarray[westwest].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());

		float vxC = Uarray[right].px()/mC;
		float vxWW = Uarray[westwest].px()/mWW;
		float vxW = Uarray[W].px()/mW;
		float vxN = Uarray[N].px()/mN;
		float vxS = Uarray[S].px()/mS;

		Uarray[right].dyvx() = (vxN - vxS) / (2.0f * dy); //OK
		Uarray[right].dxvx() = (3.0f * vxC- 4.0f * vxW + vxWW) /(2.0f * dx); //backwar difference
	}
}


void DiffOperator::VelocityXGradient_corners(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	int kp;
	float dx=fluid_class.dx,dy=fluid_class.dy;
// i=0 j=0 forward x forward y
	kp = 0 + 0 * size_x;
	float mC = fluid_class.DensityToMass(Uarray[kp].n());
	float mE = fluid_class.DensityToMass(Uarray[kp+1].n());
	float mEE = fluid_class.DensityToMass(Uarray[kp+2].n());
	float mS = fluid_class.DensityToMass(Uarray[kp+stride].n());
	float mSS = fluid_class.DensityToMass(Uarray[kp+2*stride].n());
	float vxC = Uarray[kp].px()/mC;
	float vxE = Uarray[kp+1].px()/mE;
	float vxEE = Uarray[kp+2].px()/mEE;
	float vxS = Uarray[kp+stride].px()/mS;
	float vxSS = Uarray[kp+2*stride].px()/mSS;

	Uarray[kp].dxvx() = (-3.0f * vxC + 4.0f * vxE - vxEE ) / (2.0f * dx);
	Uarray[kp].dyvx() = (-3.0f * vxC + 4.0f * vxS - vxSS ) / (2.0f * dy);
//-----------------------------------------------------------------------
// i=(size_x-1) j=0 backward x forward y
	kp = (size_x - 1) + 0 * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	float mW = fluid_class.DensityToMass(Uarray[kp-1].n());
	float mWW = fluid_class.DensityToMass(Uarray[kp-2].n());
	mS = fluid_class.DensityToMass(Uarray[kp+stride].n());
	mSS = fluid_class.DensityToMass(Uarray[kp+2*stride].n());
	vxC = Uarray[kp].px()/mC;
	float vxW = Uarray[kp-1].px()/mW;
	float vxWW = Uarray[kp-2].px()/mWW;
	vxS = Uarray[kp+stride].px()/mS;
	vxSS = Uarray[kp+2*stride].px()/mSS;

	Uarray[kp].dxvx() = (3.0f * vxC - 4.0f * vxW + vxWW ) / (2.0f * dx);
	Uarray[kp].dyvx() = (-3.0f * vxC + 4.0f * vxS - vxSS ) / (2.0f * dy);
//-----------------------------------------------------------------------

// i=0 j=(size_y-1) forward x backward y
	kp = 0 + (size_y - 1) * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	mE = fluid_class.DensityToMass(Uarray[kp+1].n());
	mEE = fluid_class.DensityToMass(Uarray[kp+2].n());
	float mN = fluid_class.DensityToMass(Uarray[kp-stride].n());
	float mNN = fluid_class.DensityToMass(Uarray[kp-2*stride].n());
	vxC = Uarray[kp].px()/mC;
	vxE = Uarray[kp+1].px()/mE;
	vxEE = Uarray[kp+2].px()/mEE;
	float vxN = Uarray[kp-stride].px()/mN;
	float vxNN = Uarray[kp-2*stride].px()/mNN;
	Uarray[kp].dxvx() = (-3.0f * vxC + 4.0f * vxE - vxEE ) / (2.0f * dx);
	Uarray[kp].dyvx() = (3.0f * vxC - 4.0f * vxN + vxNN ) / (2.0f * dy);

//-----------------------------------------------------------------------

// i=(size_x-1) j=(size_y-1) backward x backward y
	kp = (size_x - 1) + (size_y - 1) * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	mW = fluid_class.DensityToMass(Uarray[kp-1].n());
	mWW = fluid_class.DensityToMass(Uarray[kp-2].n());
	mN = fluid_class.DensityToMass(Uarray[kp-stride].n());
	mNN = fluid_class.DensityToMass(Uarray[kp-2*stride].n());

	vxC = Uarray[kp].px()/mC;
	vxW = Uarray[kp-1].px()/mW;
	vxWW = Uarray[kp-2].px()/mWW;
	vxN = Uarray[kp-stride].px()/mN;
	vxNN = Uarray[kp-2*stride].px()/mNN;

	Uarray[kp].dyvx() = (3.0f * vxC - 4.0f * vxN + vxNN ) / (2.0f * dy);
	Uarray[kp].dxvx() = (3.0f * vxC - 4.0f * vxW + vxWW ) / (2.0f * dx);
}



void DiffOperator::VelocityYGradient_bulk(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	float dx= fluid_class.dx, dy=fluid_class.dy;
#pragma omp parallel for default(none) shared(fluid_class,Uarray,stride,size_x,size_y,dx,dy)
	for (int kp = 1 + size_x; kp <= size_x * size_y - size_x - 2; kp++) {
		int N,S,E,W;
		if (kp % stride != stride - 1 && kp % stride != 0) {
			N=kp+stride;
			S=kp-stride;
			E=kp+1;
			W=kp-1;
			float mE = fluid_class.DensityToMass(Uarray[E].n());
			float mW = fluid_class.DensityToMass(Uarray[W].n());
			float mN = fluid_class.DensityToMass(Uarray[N].n());
			float mS = fluid_class.DensityToMass(Uarray[S].n());


			float vyE = Uarray[E].py()/mE;
			float vyW = Uarray[W].py()/mW;
			float vyN = Uarray[N].py()/mN;
			float vyS = Uarray[S].py()/mS;

			Uarray[kp].dxvy() = (vyE - vyW) / (2.0f * dx);
			Uarray[kp].dyvy() = (vyN - vyS) / (2.0f * dy);
		}
	}
}


void DiffOperator::VelocityYGradient_top(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int i = 1; i <= size_x - 2; i++) { // topo rede principal, ou seja j=(size_y - 1)
		int top = i + (size_y - 1) * stride;
		int southsouth = i + (size_y - 3) * stride;
		S=top-stride;
		E=top+1;
		W=top-1;
		float mC = fluid_class.DensityToMass(Uarray[top].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());
		float mSS = fluid_class.DensityToMass(Uarray[southsouth].n());

		float vyC = Uarray[top].py()/mC;
		float vyE = Uarray[E].py()/mE;
		float vyW = Uarray[W].py()/mW;
		float vySS = Uarray[southsouth].py()/mSS;
		float vyS = Uarray[S].py()/mS;

		Uarray[top].dxvy() = (vyE - vyW) / (2.0f * dx);
		Uarray[top].dyvy() = (3.0f * vyC - 4.0f * vyS + vySS) /(2.0f * dy); //backward finite difference
	}
}

void DiffOperator::VelocityYGradient_bottom(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int i = 1; i <= size_x - 2; i++) { // fundo rede principal, ou seja j=0
		int bottom = i; //i+0*nx
		int northnorth = i + 2 * stride;
		N=bottom+stride;
		E=bottom+1;
		W=bottom-1;
		float mC = fluid_class.DensityToMass(Uarray[bottom].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mNN = fluid_class.DensityToMass(Uarray[northnorth].n());


		float vyC = Uarray[bottom].py()/mC;
		float vyE = Uarray[E].py()/mE;
		float vyW = Uarray[W].py()/mW;
		float vyNN = Uarray[northnorth].py()/mNN;
		float vyN = Uarray[N].py()/mN;

		Uarray[bottom].dxvy() = (vyE - vyW) / (2.0f * dx);
		Uarray[bottom].dyvy() = (-3.0f * vyC + 4.0f * vyN - vyNN) /(2.0f * dy); //backward finite difference
	}
}

void DiffOperator::VelocityYGradient_left(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int j = 1; j <= size_y - 2; j++) { //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j * stride;
		int easteast = left + 2;
		N=left+stride;
		S=left-stride;
		E=left+1;

		float mC = fluid_class.DensityToMass(Uarray[left].n());
		float mE = fluid_class.DensityToMass(Uarray[E].n());
		float mEE = fluid_class.DensityToMass(Uarray[easteast].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());


		float vyC = Uarray[left].py()/mC;
		float vyE = Uarray[E].py()/mE;
		float vyEE = Uarray[easteast].py()/mEE;
		float vyN = Uarray[N].py()/mN;
		float vyS = Uarray[S].py()/mS;

		Uarray[left].dxvy() = (-3.0f * vyC+ 4.0f * vyE - vyEE) /(2.0f * dx); //forward difference
		Uarray[left].dyvy() = (vyN - vyS) / (2.0f * dy); //OK
	}
}

void DiffOperator::VelocityYGradient_right(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	float dx=fluid_class.dx,dy=fluid_class.dy;
	for (int j = 1; j <= size_y - 2; j++) { //lado direito da rede principal ou seja i=(size_x-1)
		int right = (size_x - 1) + j * stride;
		int westwest = right - 2;
		N=right+stride;
		S=right-stride;
		W=right-1;

		float mC = fluid_class.DensityToMass(Uarray[right].n());
		float mWW = fluid_class.DensityToMass(Uarray[westwest].n());
		float mW = fluid_class.DensityToMass(Uarray[W].n());
		float mN = fluid_class.DensityToMass(Uarray[N].n());
		float mS = fluid_class.DensityToMass(Uarray[S].n());

		float vyC = Uarray[right].py()/mC;
		float vyWW = Uarray[westwest].py()/mWW;
		float vyW = Uarray[W].py()/mW;
		float vyN = Uarray[N].py()/mN;
		float vyS = Uarray[S].py()/mS;

		Uarray[right].dyvy() = (vyN - vyS) / (2.0f * dy); //OK
		Uarray[right].dxvy() = (3.0f * vyC- 4.0f * vyW + vyWW) /(2.0f * dx);//backwar difference
	}
}

void DiffOperator::VelocityYGradient_corners(Fluid2D &fluid_class, StateVec2D *Uarray, int size_x, int size_y) {
	int stride = size_x;
	int N,S,E,W;
	int kp;
	float dx=fluid_class.dx,dy=fluid_class.dy;
// i=0 j=0 forward x forward y
	kp = 0 + 0 * size_x;
	float mC = fluid_class.DensityToMass(Uarray[kp].n());
	float mE = fluid_class.DensityToMass(Uarray[kp+1].n());
	float mEE = fluid_class.DensityToMass(Uarray[kp+2].n());
	float mS = fluid_class.DensityToMass(Uarray[kp+stride].n());
	float mSS = fluid_class.DensityToMass(Uarray[kp+2*stride].n());


	float vyC = Uarray[kp].py()/mC;
	float vyE = Uarray[kp+1].py()/mE;
	float vyEE = Uarray[kp+2].py()/mEE;
	float vyS = Uarray[kp+stride].py()/mS;
	float vySS = Uarray[kp+2*stride].py()/mSS;
	Uarray[kp].dxvy() = (-3.0f * vyC + 4.0f * vyE - vyEE ) / (2.0f * dx);
	Uarray[kp].dyvy() = (-3.0f * vyC + 4.0f * vyS - vySS ) / (2.0f * dy);
//-----------------------------------------------------------------------
// i=(size_x-1) j=0 backward x forward y
	kp = (size_x - 1) + 0 * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	float mW = fluid_class.DensityToMass(Uarray[kp-1].n());
	float mWW = fluid_class.DensityToMass(Uarray[kp-2].n());
	mS = fluid_class.DensityToMass(Uarray[kp+stride].n());
	mSS = fluid_class.DensityToMass(Uarray[kp+2*stride].n());

	vyC = Uarray[kp].py()/mC;
	float vyW = Uarray[kp-1].py()/mW;
	float vyWW = Uarray[kp-2].py()/mWW;
	vyS = Uarray[kp+stride].py()/mS;
	vySS = Uarray[kp+2*stride].py()/mSS;
	Uarray[kp].dxvy() = (3.0f * vyC - 4.0f * vyW + vyWW ) / (2.0f * dx);
	Uarray[kp].dyvy() = (-3.0f * vyC + 4.0f * vyS - vySS ) / (2.0f * dy);
//-----------------------------------------------------------------------

// i=0 j=(size_y-1) forward x backward y
	kp = 0 + (size_y - 1) * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	mE = fluid_class.DensityToMass(Uarray[kp+1].n());
	mEE = fluid_class.DensityToMass(Uarray[kp+2].n());
	float mN = fluid_class.DensityToMass(Uarray[kp-stride].n());
	float mNN = fluid_class.DensityToMass(Uarray[kp-2*stride].n());

	vyC = Uarray[kp].py()/mC;
	vyE = Uarray[kp+1].py()/mE;
	vyEE = Uarray[kp+2].py()/mEE;
	float vyN = Uarray[kp-stride].py()/mN;
	float vyNN = Uarray[kp-2*stride].py()/mNN;
	Uarray[kp].dxvy() = (-3.0f * vyC + 4.0f * vyE - vyEE ) / (2.0f *dx);
	Uarray[kp].dyvy() = (3.0f * vyC - 4.0f * vyN + vyNN ) / (2.0f * dy);

//-----------------------------------------------------------------------

// i=(size_x-1) j=(size_y-1) backward x backward y
	kp = (size_x - 1) + (size_y - 1) * size_x;
	mC = fluid_class.DensityToMass(Uarray[kp].n());
	mW = fluid_class.DensityToMass(Uarray[kp-1].n());
	mWW = fluid_class.DensityToMass(Uarray[kp-2].n());
	mN = fluid_class.DensityToMass(Uarray[kp-stride].n());
	mNN = fluid_class.DensityToMass(Uarray[kp-2*stride].n());


	vyC = Uarray[kp].py()/mC;
	vyW = Uarray[kp-1].py()/mW;
	vyWW = Uarray[kp-2].py()/mWW;
	vyN = Uarray[kp-stride].py()/mN;
	vyNN = Uarray[kp-2*stride].py()/mNN;
	Uarray[kp].dyvy() = (3.0f * vyC - 4.0f * vyN + vyNN ) / (2.0f * dy);
	Uarray[kp].dxvy() = (3.0f * vyC - 4.0f * vyW + vyWW ) / (2.0f * dx);
}

/*
DiffOperator::DiffOperator(Fluid2D &fluid) {
fluid_pointer=&fluid;
}
*/
