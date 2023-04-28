//
// Created by pcosme on 05/02/2022.
//

#include "StateVec1DLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec1D::StateVec1D(const StateVec1D & obj) {
density=obj.density;
velocity=obj.velocity;
velocity_gradient=obj.velocity_gradient;
sound=obj.sound;
}


StateVec1D::StateVec1D(float den, float vel, float snd) {
	density=den;
	velocity=vel;
	velocity_gradient=0.0f;
	sound=snd;
}


StateVec1D::StateVec1D(float den, float vel) {
	density=den;
	velocity=vel;
	velocity_gradient=0.0f;
	sound=1.0f;
}


float &StateVec1D::n() {
	return density;
}

float &StateVec1D::v() {
	return velocity;
}

float &StateVec1D::grad_v() {
	return velocity_gradient;
}

float &StateVec1D::S() {
	return sound;
}

float &StateVec1D::tmp() {
    return temperature;
}

float &StateVec1D::phi() {
	return potencial;
}



StateVec1D &StateVec1D::operator=(const StateVec1D & obj) {
	if(this != &obj) {
		this->density=obj.density;
		this->velocity=obj.velocity;
		this->sound=obj.sound;
		this->velocity_gradient=obj.velocity_gradient;
	}
	return *this;
}

StateVec1D StateVec1D::operator+(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density =this->density + obj.density;
	res.velocity = this->velocity + obj.velocity;
	res.velocity_gradient = this->velocity_gradient + obj.velocity_gradient;
	return res;
}

StateVec1D StateVec1D::operator-(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density = this->density - obj.density;
	res.velocity = this->velocity - obj.velocity;
	res.velocity_gradient = this->velocity_gradient - obj.velocity_gradient;
	return res;
}

StateVec1D StateVec1D::operator*(const StateVec1D &obj) const {
	StateVec1D res{};
	res.density = this->density * obj.density;
	res.velocity = this->velocity * obj.velocity;
	res.velocity_gradient = this->velocity_gradient * obj.velocity_gradient;
	return res;
}

StateVec1D StateVec1D::operator/(const StateVec1D &obj) const{
	StateVec1D res{};
	res.density = this->density / obj.density;
	res.velocity = this->velocity / obj.velocity;
	res.velocity_gradient = this->velocity_gradient / obj.velocity_gradient;
	return res;
}


StateVec1D operator*(const StateVec1D &obj, float value){
	StateVec1D res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	res.velocity_gradient = value*obj.velocity_gradient;
	return res;
}
StateVec1D operator*(float value, const StateVec1D &obj) {
	StateVec1D res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	res.velocity_gradient = value*obj.velocity_gradient;
	return res;
}

StateVec1D operator/(const StateVec1D &obj, float value){
	StateVec1D res{};
	res.density = obj.density/value;
	res.velocity =obj.velocity/value;
	res.velocity_gradient =obj.velocity_gradient/value;
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec1D &obj) {
	return outstream << obj.density <<"\t"<< obj.velocity <<"\t"<< obj.sound;
}



