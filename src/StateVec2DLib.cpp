//
// Created by pcosme on 05/02/2022.
//

#include "StateVec2DLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec2D::StateVec2D(const StateVec2D & obj) {
	density=obj.density;
	momentum_x=obj.momentum_x;
	momentum_y=obj.momentum_y;
	temperature=obj.temperature;
}


StateVec2D::StateVec2D(float den, float velx, float vely,float temp) {
	density=den;
	momentum_x=velx;
	momentum_y=vely;
	temperature=temp;
}

StateVec2D::StateVec2D(float den, float velx, float vely,float temp,float snd) {
	density=den;
	momentum_x=velx;
	momentum_y=vely;
	temperature=temp;
	sound=snd;
}


float &StateVec2D::n() {
	return density;
}

float &StateVec2D::px() {
	return momentum_x;
}

float &StateVec2D::py() {
	return momentum_y;
}

float &StateVec2D::tmp() {
	return temperature;
}

float &StateVec2D::S() {
	return sound;
}


StateVec2D &StateVec2D::operator=(const StateVec2D & obj) {
	if(this != &obj) {
		density=obj.density;
		momentum_x=obj.momentum_x;
		momentum_y=obj.momentum_y;
		temperature=obj.temperature;
		sound=obj.sound;
	}
	return *this;
}

StateVec2D StateVec2D::operator+(const StateVec2D &obj) const {
	StateVec2D res{};
	res.density = density + obj.density;
	res.momentum_x = momentum_x + obj.momentum_x;
	res.momentum_y = momentum_y + obj.momentum_y;
	res.temperature = temperature + obj.temperature;
	res.sound = sound + obj.sound;
	return res;
}

StateVec2D StateVec2D::operator-(const StateVec2D &obj) const {
	StateVec2D res{};
	res.density = density - obj.density;
	res.momentum_x = momentum_x - obj.momentum_x;
	res.momentum_y = momentum_y - obj.momentum_y;
	res.temperature = temperature - obj.temperature;
	res.sound = sound - obj.sound;
	return res;
}

StateVec2D StateVec2D::operator*(const StateVec2D &obj) const {
	StateVec2D res{};
	res.density = density * obj.density;
	res.momentum_x = momentum_x * obj.momentum_x;
	res.momentum_y = momentum_y * obj.momentum_y;
	res.temperature = temperature * obj.temperature;
	res.sound = sound * obj.sound;
	return res;
}

StateVec2D StateVec2D::operator/(const StateVec2D &obj) const{
	StateVec2D res{};
	res.density = density / obj.density;
	res.momentum_x = momentum_x / obj.momentum_x;
	res.momentum_y = momentum_y / obj.momentum_y;
	res.temperature = temperature / obj.temperature;
	res.sound = sound / obj.sound;
	return res;
}


StateVec2D operator*(const StateVec2D &obj, float value){
	StateVec2D res{};
	res.density = value*obj.density;
	res.momentum_x = value * obj.momentum_x;
	res.momentum_y = value * obj.momentum_y;
	res.temperature = value*obj.temperature;
	res.sound = value*obj.sound;
	return res;
}
StateVec2D operator*(float value, const StateVec2D &obj) {
	StateVec2D res{};
	res.density = value*obj.density;
	res.momentum_x = value * obj.momentum_x;
	res.momentum_y = value * obj.momentum_y;
	res.temperature = value*obj.temperature;
	res.sound = value*obj.sound;
	return res;
}

StateVec2D operator/(const StateVec2D &obj, float value){
	StateVec2D res{};
	res.density = obj.density/value;
	res.momentum_x = obj.momentum_x / value;
	res.momentum_y = obj.momentum_y / value;
	res.temperature =obj.temperature/value;
	res.sound =obj.sound/value;
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec2D &obj) {
	return outstream << obj.density << "\t" << obj.momentum_x << "\t" << obj.momentum_y << "\t" << obj.temperature << "\t" << obj.sound;
}


