//
// Created by pcosme on 27/12/2020.
//

#ifndef RICHTMYER2DLIB_H
#define RICHTMYER2DLIB_H


#include "TethysBaseLib.h"

class Richtmyer2D{
public:
	float * Data;

	void RichtmyerFirstStep();
	void RichtmyerSecondStep();

	virtual float FluxX(float n, float flx_x, float flx_y, float mass, float s ) = 0;
	virtual float FluxY(float n, float flx_x, float flx_y, float mass, float s ) = 0;
	virtual float Source(float n, float flx_x, float flx_y, float mass, float s) = 0;

};


#endif //RICHTMYER2DLIB_H

