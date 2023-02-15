/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef STATEVEC1DLIB_H
#define STATEVEC1DLIB_H

#include <iostream>

class StateVec1D{
private:
	float density;
	float velocity;
public:
	StateVec1D()=default; //default constructor
	StateVec1D(const StateVec1D&); //copy constructor
	StateVec1D(float den, float vel);
	~StateVec1D()=default;

	float& v();
	float& n();

	StateVec1D& operator=(const StateVec1D&);
	StateVec1D operator + (StateVec1D const &obj) const ;
	StateVec1D operator - (StateVec1D const &obj) const ;
	StateVec1D operator * (StateVec1D const &obj) const ;
	StateVec1D operator / (StateVec1D const &obj) const ;
	friend StateVec1D operator*(const StateVec1D &obj, float value);
	friend StateVec1D operator*(float value, const StateVec1D &obj);
	friend StateVec1D operator/(const StateVec1D &obj, float value);
	friend std::ostream& operator<<(std::ostream& outstream, const StateVec1D &obj);
};




#endif //STATEVEC1DLIB_H


