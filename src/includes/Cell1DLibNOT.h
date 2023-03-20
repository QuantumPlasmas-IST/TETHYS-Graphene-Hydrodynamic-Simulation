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
	StateVec1D *U_ptr;
	Fluid1D *fluid_ptr;
public:
	CellHandler1D(int, int, Fluid1D *, StateVec1D *);
	CellHandler1D(int, int, StateVec1D *);
	~CellHandler1D()=default;

//	StateVec1D TVD(char side,char edge);


//TODO adicionar mais flux limiters

//	StateVec1D VanLeer(int i);
//	StateVec1D Roe(int i);

};



class NumericalFlux{
public:
	static StateVec1D Average(Fluid1D* fluido, StateVec1D L, StateVec1D R);
//	static StateVec1D Central(Fluid1D* fluido, StateVec1D L, StateVec1D R);
//	static StateVec1D Characteristic(Fluid1D* fluido, StateVec1D L, StateVec1D R);

};

class NumericalSource{
public:
	static StateVec1D Average(Fluid1D* fluido, StateVec1D L, StateVec1D R);
};

#endif //CELL1D_H
