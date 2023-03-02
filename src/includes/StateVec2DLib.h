/************************************************************************************************\
* 2022 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* 																 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef STATEVECLIBFERMI2D_H
#define STATEVECLIBFERMI2D_H

#include <iostream>

class StateVec2D{
private:
	float density;
	float momentum_x;
	float momentum_y;
	float sound=1.0f;
	float temperature;
	float velXGradient_x;
	float velXGradient_y;
	float velYGradient_x;
	float velYGradient_y;
	float velXLaplacian;
	float velYLaplacian;
	float tmpLaplacian;
public:
	StateVec2D()=default; //default constructor
	StateVec2D(const StateVec2D&); //copy constructor
	StateVec2D(float den, float velx, float vely, float temp);
	StateVec2D(float den, float velx, float vely, float temp, float snd);
	~StateVec2D()=default;

	float& px();
	float& py();
	float& n();
	float& tmp();
	float& S();
	float& d2vx();
	float& d2vy();
	float& d2tmp();

	float& dxvx();
	float& dyvx();
	float& dxvy();
	float& dyvy();

	StateVec2D& operator=(const StateVec2D&);
	StateVec2D operator + (StateVec2D const &obj) const ;
	StateVec2D operator - (StateVec2D const &obj) const ;
	StateVec2D operator * (StateVec2D const &obj) const ;
	StateVec2D operator / (StateVec2D const &obj) const ;
	friend StateVec2D operator*(const StateVec2D &obj, float value);
	friend StateVec2D operator*(float value, const StateVec2D &obj);
	friend StateVec2D operator/(const StateVec2D &obj, float value);
	friend std::ostream& operator<<(std::ostream& outstream, const StateVec2D &obj);
};




#endif //STATEVECLIBFERMI2D_H


