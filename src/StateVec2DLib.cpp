//
// Created by pcosme on 05/02/2022.
//

#include "StateVec2DLib.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StateVec2D::StateVec2D(const StateVec2D & obj) {
	density=obj.density;
	momentum_x=obj.momentum_x;
	momentum_y=obj.momentum_y;

	//velXGradient_x = obj.velXGradient_x ;
	//velXGradient_y = obj.velXGradient_y;
	//velYGradient_x = obj.velYGradient_x;
	//velYGradient_y = obj.velYGradient_y;s
	velXLaplacian = obj.velXLaplacian;
	velYLaplacian = obj.velYLaplacian;
	tmpLaplacian = obj.tmpLaplacian;

	temperature=obj.temperature;
	sound=obj.sound;
}



StateVec2D::StateVec2D(float den, float velx, float vely) {
	density=den;
	momentum_x=velx;
	momentum_y=vely;
	velXGradient_x = 0.0f;
	velXGradient_y = 0.0f;
	velYGradient_x = 0.0f;
	velYGradient_y = 0.0f;
	velXLaplacian = 0.0f;
	velYLaplacian = 0.0f;
	tmpLaplacian = 0.0f;
}

StateVec2D::StateVec2D(float den, float velx, float vely, float temp) : StateVec2D(den, velx, vely)  {
	this->temperature=temp;
}

StateVec2D::StateVec2D(float den, float velx, float vely, float temp, float snd) : StateVec2D(den, velx, vely, temp) {
	this->sound=snd;
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

float &StateVec2D::d2vx(){
	return velXLaplacian;
}

float &StateVec2D::d2vy(){
	return velYLaplacian;
}

float &StateVec2D::d2tmp(){
	return tmpLaplacian;
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

		velXGradient_x=obj.velXGradient_x;
		velXGradient_y=obj.velXGradient_y;
		velYGradient_x=obj.velYGradient_x;
		velYGradient_y=obj.velYGradient_y;

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

	res.velXGradient_x=velXGradient_x+obj.velXGradient_x;
	res.velXGradient_y=velXGradient_y+obj.velXGradient_y;
	res.velYGradient_x=velYGradient_x+obj.velYGradient_x;
	res.velYGradient_y=velYGradient_y+obj.velYGradient_y;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}

StateVec2D StateVec2D::operator-(const StateVec2D &obj) const {
	StateVec2D res{};
	res.density = density - obj.density;
	res.momentum_x = momentum_x - obj.momentum_x;
	res.momentum_y = momentum_y - obj.momentum_y;
	res.temperature = temperature - obj.temperature;

	res.velXGradient_x=velXGradient_x-obj.velXGradient_x;
	res.velXGradient_y=velXGradient_y-obj.velXGradient_y;
	res.velYGradient_x=velYGradient_x-obj.velYGradient_x;
	res.velYGradient_y=velYGradient_y-obj.velYGradient_y;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}

StateVec2D StateVec2D::operator*(const StateVec2D &obj) const {
	StateVec2D res{};
	res.density = density * obj.density;
	res.momentum_x = momentum_x * obj.momentum_x;
	res.momentum_y = momentum_y * obj.momentum_y;
	res.temperature = temperature * obj.temperature;

	res.velXGradient_x=velXGradient_x*obj.velXGradient_x;
	res.velXGradient_y=velXGradient_y*obj.velXGradient_y;
	res.velYGradient_x=velYGradient_x*obj.velYGradient_x;
	res.velYGradient_y=velYGradient_y*obj.velYGradient_y;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}

StateVec2D StateVec2D::operator/(const StateVec2D &obj) const{
	StateVec2D res{};
	res.density = density / obj.density;
	res.momentum_x = momentum_x / obj.momentum_x;
	res.momentum_y = momentum_y / obj.momentum_y;
	res.temperature = temperature / obj.temperature;


	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}


StateVec2D operator*(const StateVec2D &obj, float value){
	StateVec2D res{};
	res.density = value*obj.density;
	res.momentum_x = value * obj.momentum_x;
	res.momentum_y = value * obj.momentum_y;
	res.temperature = value*obj.temperature;

	res.velXGradient_x=value*obj.velXGradient_x;
	res.velXGradient_y=value*obj.velXGradient_y;
	res.velYGradient_x=value*obj.velYGradient_x;
	res.velYGradient_y=value*obj.velYGradient_y;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}
StateVec2D operator*(float value, const StateVec2D &obj) {
	StateVec2D res{};
	res.density = value*obj.density;
	res.momentum_x = value * obj.momentum_x;
	res.momentum_y = value * obj.momentum_y;
	res.temperature = value*obj.temperature;

	res.velXGradient_x=value*obj.velXGradient_x;
	res.velXGradient_y=value*obj.velXGradient_y;
	res.velYGradient_x=value*obj.velYGradient_x;
	res.velYGradient_y=value*obj.velYGradient_y;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}

StateVec2D operator/(const StateVec2D &obj, float value){
	StateVec2D res{};
	res.density = obj.density/value;
	res.momentum_x = obj.momentum_x / value;
	res.momentum_y = obj.momentum_y / value;
	res.temperature =obj.temperature/value;

	res.velXGradient_x=obj.velXGradient_x / value;
	res.velXGradient_y=obj.velXGradient_y / value;
	res.velYGradient_x=obj.velYGradient_x / value;
	res.velYGradient_y=obj.velYGradient_y / value;

	res.sound = obj.sound; //simplesmente atribuir o valore de vel do som do obj
	return res;
}

std::ostream &operator<<(std::ostream &outstream, const StateVec2D &obj) {
	return outstream << obj.density << "\t" << obj.momentum_x << "\t" << obj.momentum_y << "\t" << obj.temperature << "\t" << obj.sound;
}

float &StateVec2D::dxvx() {
	return velXGradient_x;
}

float &StateVec2D::dyvx() {
	return velXGradient_y;
}

float &StateVec2D::dxvy() {
	return velYGradient_x;
}

float &StateVec2D::dyvy() {
	return velYGradient_y;
}


