/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef CELL2D_H
#define CELL2D_H

#include "includes/Fluid2DLib.h"



class CellHandler2D {
private:
	int index;
	StateVec2D U{};
	StateVec2D *U_ptr;
	Fluid2D *fluid_ptr;
public:
	CellHandler2D(int, Fluid2D *, StateVec2D *);
	~CellHandler2D()=default;
//	StateVec &E(StateVec * Uin);
//	StateVec &W(StateVec * Uin);

//	StateVec2D TVD(char side,char edge);
//	StateVec2D UNO(char side,char edge); //TODO implementar a reconstrução com o UNO2

//	float VanLeer(int i);
//	float Roe(int i); //TODO adicionar mais flux limiters
};



class NumericalFlux{
public:
	static StateVec2D Average(Fluid2D* fluido, StateVec2D L, StateVec2D R);
	static StateVec2D Central(Fluid2D* fluido, StateVec2D L, StateVec2D R);
	static StateVec2D Characteristic(Fluid2D* fluido, StateVec2D L, StateVec2D R);

};

#endif CELL2D_H
