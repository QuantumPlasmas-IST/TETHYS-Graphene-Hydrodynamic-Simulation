/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef STATEVECLIB_H
#define STATEVECLIB_H

#include <iostream>

class StateVec{
private:
	float density;
	float velocity;
	float sound=1.0f;
public:
	StateVec()=default; //default constructor
	StateVec(const StateVec&); //copy constructor
	StateVec(float den,float vel,float snd);
	StateVec(float den,float vel);
	~StateVec()=default;



	float& v();
	float& n();
	float& S();

	StateVec& operator=(const StateVec&);
	StateVec operator + (StateVec const &obj) const ;
	StateVec operator - (StateVec const &obj) const ;
	StateVec operator * (StateVec const &obj) const ;
	StateVec operator / (StateVec const &obj) const ;
	friend StateVec operator*(const StateVec &obj, float value);
	friend StateVec operator*(float value, const StateVec &obj);
	friend StateVec operator/(const StateVec &obj, float value);
	friend std::ostream& operator<<(std::ostream& outstream, const StateVec &obj);
};

#endif //STATEVECLIB_H


