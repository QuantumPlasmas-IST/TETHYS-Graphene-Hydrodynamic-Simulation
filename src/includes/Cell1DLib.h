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
	StateVec1D U{};
	StateVec1D *U_ptr;
	Fluid1D *fluid_ptr;
public:
	CellHandler1D(int, Fluid1D *, StateVec1D *);
	~CellHandler1D()=default;
//	StateVec1D &E(StateVec1D * Uin);
//	StateVec1D &W(StateVec1D * Uin);

	StateVec1D TVD(char side, char edge);
	StateVec1D UNO(char side, char edge); //TODO implementar a reconstrução com o UNO2

	float VanLeer(int i);
	float Roe(int i); //TODO adicionar mais flux limiters
};



class NumericalFlux{
public:
	static StateVec1D Average(Fluid1D* fluido, StateVec1D L, StateVec1D R);
	static StateVec1D Central(Fluid1D* fluido, StateVec1D L, StateVec1D R);
	static StateVec1D Characteristic(Fluid1D* fluido, StateVec1D L, StateVec1D R);

};

#endif //CELL1D_H
