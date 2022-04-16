//
// Created by pcosme on 05/02/2022.
//

#include "StateVecLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec::StateVec(const StateVec & obj) {
density=obj.density;
density_laplacian=obj.density_laplacian;
velocity=obj.velocity;
sound=obj.sound;
}


StateVec::StateVec(float den, float vel,float snd) {
	density=den;
	density_laplacian=0.0f;
	velocity=vel;
	sound=snd;
}


StateVec::StateVec(float den, float vel) {
	density=den;
	density_laplacian=0.0f;
	velocity=vel;
	sound=1.0f;
}


float &StateVec::n() {
	return density;
}

float &StateVec::lap_n() {
	return density_laplacian;
}

float &StateVec::v() {
	return velocity;
}

float &StateVec::S() {
	return sound;
}


StateVec &StateVec::operator=(const StateVec & obj) {
	if(this != &obj) {
		this->density=obj.density;
		this->density_laplacian=obj.density_laplacian;
		this->velocity=obj.velocity;
		this->sound=obj.sound;
	}
	return *this;
}

StateVec StateVec::operator+(const StateVec &obj) const {
	StateVec res{};
	res.density =this->density + obj.density;
	res.density_laplacian = this->density_laplacian + obj.density_laplacian;
	res.velocity = this->velocity + obj.velocity;
	return res;
}

StateVec StateVec::operator-(const StateVec &obj) const {
	StateVec res{};
	res.density = this->density - obj.density;
	res.density_laplacian = this->density_laplacian - obj.density_laplacian;
	res.velocity = this->velocity - obj.velocity;
	return res;
}

StateVec StateVec::operator*(const StateVec &obj) const {
	StateVec res{};
	res.density = this->density * obj.density;
	res.density_laplacian = this->density_laplacian * obj.density_laplacian;
	res.velocity = this->velocity * obj.velocity;
	return res;
}

StateVec StateVec::operator/(const StateVec &obj) const{
	StateVec res{};
	res.density = this->density / obj.density;
	res.density_laplacian = this->density_laplacian / obj.density_laplacian;
	res.velocity = this->velocity / obj.velocity;
	return res;
}


StateVec operator*(const StateVec &obj, float value){
	StateVec res{};
	res.density = value*obj.density;
	res.density_laplacian = value*obj.density_laplacian;
	res.velocity = value*obj.velocity;
	return res;
}
StateVec operator*(float value, const StateVec &obj) {
	StateVec res{};
	res.density = value*obj.density;
	res.density_laplacian = value*obj.density_laplacian;
	res.velocity = value*obj.velocity;
	return res;
}

StateVec operator/(const StateVec &obj, float value){
	StateVec res{};
	res.density = obj.density/value;
	res.density_laplacian = obj.density_laplacian/value;
	res.velocity =obj.velocity/value;
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec &obj) {
	return outstream << obj.density <<"\t"<< obj.velocity <<"\t"<< obj.sound;
}





