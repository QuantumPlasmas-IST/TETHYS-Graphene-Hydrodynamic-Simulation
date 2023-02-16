//
// Created by pcosme on 05/02/2022.
//

#include "StateVecLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec::StateVec(const StateVec & obj) {
density=obj.density;
velocity=obj.velocity;
velocity_gradient=obj.velocity_gradient;
sound=obj.sound;
}


StateVec::StateVec(float den, float vel,float snd) {
	density=den;
	velocity=vel;
	velocity_gradient=0.0f;
	sound=snd;
}


StateVec::StateVec(float den, float vel) {
	density=den;
	velocity=vel;
	velocity_gradient=0.0f;
	sound=1.0f;
}


float &StateVec::n() {
	return density;
}

float &StateVec::v() {
	return velocity;
}

float &StateVec::grad_v() {
	return velocity_gradient;
}

float &StateVec::S() {
	return sound;
}

StateVec &StateVec::operator=(const StateVec & obj) {
	if(this != &obj) {
		this->density=obj.density;
		this->velocity=obj.velocity;
		this->sound=obj.sound;
		this->velocity_gradient=obj.velocity_gradient;
	}
	return *this;
}

StateVec StateVec::operator+(const StateVec &obj) const {
	StateVec res{};
	res.density =this->density + obj.density;
	res.velocity = this->velocity + obj.velocity;
	res.velocity_gradient = this->velocity_gradient + obj.velocity_gradient;
	return res;
}

StateVec StateVec::operator-(const StateVec &obj) const {
	StateVec res{};
	res.density = this->density - obj.density;
	res.velocity = this->velocity - obj.velocity;
	res.velocity_gradient = this->velocity_gradient - obj.velocity_gradient;
	return res;
}

StateVec StateVec::operator*(const StateVec &obj) const {
	StateVec res{};
	res.density = this->density * obj.density;
	res.velocity = this->velocity * obj.velocity;
	res.velocity_gradient = this->velocity_gradient * obj.velocity_gradient;
	return res;
}

StateVec StateVec::operator/(const StateVec &obj) const{
	StateVec res{};
	res.density = this->density / obj.density;
	res.velocity = this->velocity / obj.velocity;
	res.velocity_gradient = this->velocity_gradient / obj.velocity_gradient;
	return res;
}


StateVec operator*(const StateVec &obj, float value){
	StateVec res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	res.velocity_gradient = value*obj.velocity_gradient;
	return res;
}
StateVec operator*(float value, const StateVec &obj) {
	StateVec res{};
	res.density = value*obj.density;
	res.velocity = value*obj.velocity;
	res.velocity_gradient = value*obj.velocity_gradient;
	return res;
}

StateVec operator/(const StateVec &obj, float value){
	StateVec res{};
	res.density = obj.density/value;
	res.velocity =obj.velocity/value;
	res.velocity_gradient =obj.velocity_gradient/value;
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec &obj) {
	return outstream << obj.density <<"\t"<< obj.velocity <<"\t"<< obj.sound;
}



