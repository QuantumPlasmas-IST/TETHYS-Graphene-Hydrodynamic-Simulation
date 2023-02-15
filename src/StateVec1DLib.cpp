//
// Created by pcosme on 05/02/2022.
//

#include "StateVec1DLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec1D::StateVec1D(const StateVec1D & obj) {
density=obj.density;
velocity=obj.velocity;
}


StateVec1D::StateVec1D(float den, float vel) {
	density=den;
	velocity=vel;
}


float &StateVec1D::n() {
	return density;
}

float &StateVec1D::v() {
	return velocity;
}



StateVec1D &StateVec1D::operator=(const StateVec1D & obj) {
	if(this != &obj) {
		density=obj.density;
		velocity=obj.velocity;
	}
	return *this;
}

StateVec1D StateVec1D::operator+(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density = density + obj.density;
	res.velocity = velocity + obj.velocity;
	return res;
}

StateVec1D StateVec1D::operator-(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density = density - obj.density;
	res.velocity = velocity - obj.velocity;
	return res;
}

StateVec1D StateVec1D::operator*(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density = density * obj.density;
	res.velocity = velocity * obj.velocity;
	return res;
}

StateVec1D StateVec1D::operator/(const StateVec1D &obj) const{
	StateVec1D res{};
	res.density = density / obj.density;
	res.velocity = velocity / obj.velocity;
	return res;
}


StateVec1D operator*(const StateVec1D &obj, float value){
	StateVec1D res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	return res;
}
StateVec1D operator*(float value, const StateVec1D &obj) {
	StateVec1D res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	return res;
}

StateVec1D operator/(const StateVec1D &obj, float value){
	StateVec1D res{};
	res.density = obj.density/value;
	res.velocity =obj.velocity/value;
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec1D &obj) {
	return outstream << obj.density <<"\t"<< obj.velocity;
}


