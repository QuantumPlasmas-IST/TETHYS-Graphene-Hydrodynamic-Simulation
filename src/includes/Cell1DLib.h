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
	StateVec U{};
	StateVec *U_ptr;
	Fluid1D *fluid_ptr;
public:
	CellHandler1D(int, Fluid1D *, StateVec *);
	CellHandler1D(int, StateVec *);
	~CellHandler1D()=default;
//	StateVec &E(StateVec * Uin);
//	StateVec &W(StateVec * Uin);

	StateVec TVD(char side,char edge);
	StateVec TVD(StateVec *,int pos,char side,char edge);

	StateVec UNO(char side,char edge); //TODO implementar a reconstrução com o UNO2


	StateVec VanLeerU(int i);
	float VanLeer(int i);
	float VanLeer(StateVec*,int i);
	float Roe(int i); //TODO adicionar mais flux limiters
};



class NumericalFlux{
public:
	static StateVec Average(Fluid1D* fluido, StateVec L, StateVec R);
	static StateVec Central(Fluid1D* fluido, StateVec L, StateVec R);
	static StateVec Characteristic(Fluid1D* fluido, StateVec L, StateVec R);

};

#endif //CELL1D_H
