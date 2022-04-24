/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef CELL1D_H
#define CELL1D_H

#include "includes/Fluid1DLib.h"



class CellHandler1D {
private:
	int index;
	int size;
	StateVec U{};
	StateVec *U_ptr;
	Fluid1D *fluid_ptr;
public:
	CellHandler1D(int, int, Fluid1D *, StateVec *);
	CellHandler1D(int, int, StateVec *);
	~CellHandler1D()=default;

	StateVec TVD(char side,char edge);
	StateVec TVD(StateVec *,int pos,char side,char edge);


//TODO adicionar mais flux limiters
	StateVec VanLeer(StateVec*Uin, int i);
	StateVec VanLeer(int i);

	StateVec Roe(StateVec*Uin, int i);
	StateVec Roe(int i);

};



class NumericalFlux{
public:
	static StateVec Average(Fluid1D* fluido, StateVec L, StateVec R);
	static StateVec Central(Fluid1D* fluido, StateVec L, StateVec R);
	static StateVec Characteristic(Fluid1D* fluido, StateVec L, StateVec R);

};

#endif //CELL1D_H
