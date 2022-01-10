/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/TethysMathLib.h"
#include "TethysBaseLib.h"

using namespace std;

float MathUtils::Sound_Velocity_Anisotropy(float x, float s) {
	return s;
}
float  MathUtils::Sound_Velocity_Anisotropy(float x, float y, float s) {
	float s_mod;
	s_mod=s;
	return s_mod;
}

float  MathUtils::Integral_1_D(int n, float ds, const float * f){
	float itg=0.0;
	for(int j=1; j < n / 2; j++){
		itg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	itg = itg*ds/3.0f;
	return itg;
}

float  MathUtils::Integral_2_D(int n, int m, float dx, float dy, const float * f){
	float itg;
	float interior=0.0f;
	float edges=0.0f;
	float vertices;

	vertices = f[0] + f[n - 1] + f[n * m - n] + f[n * m - 1];

	for(int i=1; i <= n - 2; i++){
		edges += f[i] + f[i+ (m - 1) * n];
	}
	for(int j=1; j <= m - 2; j++){
		edges += f[j * n] + f[n - 1 + j * n];
	}
	for(int k= 1 + n; k <= n * m - n - 2; k++) {
		if (k % n != n - 1 && k % n != 0) {
			interior += f[k];
		}
	}

	itg = dx * dy * (0.25f * vertices + 0.5f * edges + interior);
	return itg;
}


void  MathUtils::Average_Filter(const float * vec_in, float * vec_out, int size , int width ){
	for ( int i = 0; i < size; i++ ){
		if(i>=width &&i<=size-1-width){
			for(int k = i-width; k <= i+width;k++){
				vec_out[i] += vec_in[k];
			}
			vec_out[i] = vec_out[i]/(2.0f*(float)width+1.0f);
		}
		else{
			vec_out[i] = vec_in[i] ;
		}
	}
}



float  MathUtils::Signal_Average(int n, float dt, const float * f){
	float avg=0.0;
	for(int j=1; j < n / 2; j++){
		avg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	avg = avg*dt/3.0f;
	avg = avg/(static_cast<float>(n) * dt);
	return avg;
}
float  MathUtils::Gauss_Kernel(int position , float t){
	auto pos=static_cast<float>(position);
	return exp(-0.5f * pos * pos / t) / (sqrt(2.0f * MAT_PI * t));
}

float  MathUtils::Gauss_Kernel_Derivative(int position , float t){
	auto pos=static_cast<float>(position);
	return (-pos * exp(-0.5f * pos * pos / t) / t) / (sqrt(2.0f * MAT_PI * t));
}

void  MathUtils::Convolve_Gauss(unsigned int type, unsigned int m, float t, const float * in, float * out, unsigned long size){
	if(type==0){
		for(unsigned int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-1*static_cast<int>(m); k <= static_cast<int>(m); k++){
					out[i] += in[i-k] * Gauss_Kernel(k, t);
				}
			}
		}
	}
	if(type==1){
		for(unsigned int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-1*static_cast<int>(m); k <= static_cast<int>(m); k++){
					out[i] += in[i-k] * Gauss_Kernel_Derivative(k, t);
				}
			}
		}
	}
}

float  MathUtils::Stair_Case_Function(float x, float step_width, float smoothness) {
	return (0.5f*tanh(smoothness*((-(x - 1.0f)/step_width - floor(-(x - 1.0f)/step_width)) - 0.5f))/tanh(0.5f*smoothness) + 0.5f +floor(-(x - 1.0f)/step_width));
}

void MathUtils::LaplacianField(const float *array_in, float *array_out, float dx, int size_x, int size_y) {
	int stride = size_x;
#pragma omp parallel for default(none) shared(size_x,size_y,dx,stride,array_in,array_out)
	for(int kp=1+size_x; kp<=size_x*size_y-size_x-2; kp++){
		if( kp%stride!=stride-1 && kp%stride!=0){
			GridPoint point(kp,size_x,size_y,false);
			array_out[kp] = (-4.0f*array_in[point.C] +  array_in[point.N] + array_in[point.S] + array_in[point.E] + array_in[point.W] )/(dx*dx);
		}
	}
// -2	5	-4	1
	for(int i=1 ; i<=size_x-2; i++){ // topo rede principal, ou seja j=(size_y - 1)
		int top= i + (size_y - 1) * stride;
		GridPoint point(top,size_x,size_y,false);
		int southsouth= i + (size_y - 3) * stride;
		int southsouthsouth= i + (size_y - 4) * stride;
		float aux1 = -2.0f*array_in[point.C] + array_in[point.E] + array_in[point.W]; //OK
		float aux2 = 2.0f*array_in[point.C] -5.0f*array_in[point.S] +4.0f*array_in[southsouth] -1.0f*array_in[southsouthsouth];
		array_out[top] = (aux1+aux2)/(dx*dx);
	}
// 2	−5	4	−1
	for(int i=1 ; i<=size_x-2; i++){ // fundo rede principal, ou seja j=0
		int bottom=i; //i+0*nx
		GridPoint point(bottom,size_x,size_y,false);
		int northnorth=i+2*stride;
		int northnorthnorth=i+3*stride;
		float aux1 = -2.0f*array_in[point.C] + array_in[point.E] + array_in[point.W]; //OK
		float aux2 =  2.0f*array_in[point.C] -5.0f*array_in[point.N] +4.0f*array_in[northnorth] -1.0f*array_in[northnorthnorth];
		array_out[bottom] = (aux1+aux2)/(dx*dx);
	}
// 2	−5	4	−1
	for(int j=1; j<=size_y-2;j++){ //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j*stride;
		GridPoint point(left,size_x,size_y,false);
		int easteast = left + 2;
		int easteasteast = left + 3;
		float aux1 =  -2.0f*array_in[point.C] + array_in[point.N] + array_in[point.S];
		float aux2 =  2.0f*array_in[point.C] -5.0f*array_in[point.E] +4.0f*array_in[easteast] -1.0f*array_in[easteasteast];
		array_out[left] = (aux1+aux2)/(dx*dx);
	}
// -2	5	-4	1
	for(int j=1; j<=size_y-2;j++){ //lado direito da rede principal ou seja i=(Nx-1)
		int right = (size_x-1) + j*stride;
		GridPoint point(right,size_x,size_y,false);
		int westwest = right-2;
		int westwestwest = right-3;
		float aux1 =  -2.0f*array_in[right] + array_in[right+stride] + array_in[right-stride];
		float aux2 =  2.0f*array_in[right] -5.0f*array_in[right-1] +4.0f*array_in[westwest] -1.0f*array_in[westwestwest];
		array_out[right] =   (aux1+aux2)/(dx*dx);
	}
//2	−5	4	−1
	int kp;
// i=0 j=0 forward x forward y
	kp = 0 + 0*size_x;
	float aux1 = 2.0f*array_in[kp] -5.0f*array_in[kp+1] +4.0f*array_in[kp+2] -1.0f*array_in[kp+3];
	float aux2 = 2.0f*array_in[kp] -5.0f*array_in[kp+size_x] +4.0f*array_in[kp+2*size_x] -1.0f*array_in[kp+3*size_x];
	array_out[kp] = (aux1+aux2)/(dx*dx);
// i=(Nx-1) j=0 backward x forward y
	kp = (size_x-1) + 0*size_x;
	aux1 = -2.0f*array_in[kp] +5.0f*array_in[kp-1] -4.0f*array_in[kp-2] +1.0f*array_in[kp-3];
	aux2 = 2.0f*array_in[kp] -5.0f*array_in[kp+size_x] +4.0f*array_in[kp+2*size_x] -1.0f*array_in[kp+3*size_x];
	array_out[kp] =(-aux1+aux2)/(dx*dx);
// i=0 j=(size_y-1) forward x backward y
	kp = 0 + (size_y-1)*size_x;
	aux1 = 2.0f*array_in[kp] -5.0f*array_in[kp+1] +4.0f*array_in[kp+2] -1.0f*array_in[kp+3];
	aux2 = -2.0f*array_in[kp] +5.0f*array_in[kp-size_x] -4.0f*array_in[kp-2*size_x] +1.0f*array_in[kp-3*size_x];
	array_out[kp] = (aux1-aux2)/(dx*dx);
// i=(Nx-1) j=(size_y-1) backward x backward y
	kp = (size_x-1) + (size_y-1)*size_x;
	aux1 = -2.0f*array_in[kp] +5.0f*array_in[kp-1] -4.0f*array_in[kp-2] +1.0f*array_in[kp-3];
	aux2 = -2.0f*array_in[kp] +5.0f*array_in[kp-size_x] -4.0f*array_in[kp-2*size_x] +1.0f*array_in[kp-3*size_x];
	array_out[kp] = (-aux1-aux2)/(dx*dx);
}

void
MathUtils::GradientField(const float *array_in, float *array_out_x, float *array_out_y, float dx, float dy, int size_x, int size_y) {
	int stride = size_x;
#pragma omp parallel for default(none) shared(size_x, size_y, dx, dy, stride, array_in, array_out_x, array_out_y)
	for (int kp = 1 + size_x; kp <= size_x * size_y - size_x - 2; kp++) {
		if (kp % stride != stride - 1 && kp % stride != 0) {
			GridPoint point(kp, size_x, size_y, false);
			array_out_x[kp] = (array_in[point.E] - array_in[point.W]) / (2.0f * dx);
			array_out_y[kp] = (array_in[point.N] - array_in[point.S]) / (2.0f * dy);
		}
	}

	for (int i = 1; i <= size_x - 2; i++) { // topo rede principal, ou seja j=(size_y - 1)
		int top = i + (size_y - 1) * stride;
		GridPoint point(top, size_x, size_y, false);
		int southsouth = i + (size_y - 3) * stride;
		array_out_x[top] = ((array_in[point.E]) - (array_in[point.W])) / (2.0f * dx); //OK
		array_out_y[top] = (3.0f * (array_in[point.C]) - 4.0f * (array_in[point.S]) + 1.0f * (array_in[southsouth])) /
		                   (2.0f * dy); //backward finite difference
	}
	for (int i = 1; i <= size_x - 2; i++) { // fundo rede principal, ou seja j=0
		int bottom = i; //i+0*nx
		GridPoint point(bottom, size_x, size_y, false);
		int northnorth = i + 2 * stride;
		array_out_x[bottom] = ((array_in[point.E]) - (array_in[point.W])) / (2.0f * dx);
		array_out_y[bottom] = (-3.0f * (array_in[point.C]) + 4.0f * (array_in[point.N]) - 1.0f * (array_in[northnorth])) /
		                      (2.0f * dy); //forward finite difference
	}
	for (int j = 1; j <= size_y - 2; j++) { //lado esquerdo da rede principal ou seja i=0
		int left = 0 + j * stride;
		GridPoint point(left, size_x, size_y, false);
		int easteast = left + 2;
		array_out_x[left] = (-3.0f * (array_in[point.C]) + 4.0f * (array_in[point.E]) - 1.0f * (array_in[easteast])) /
		                    (2.0f * dx); //forward difference
		array_out_y[left] = ((array_in[point.N]) - (array_in[point.S])) / (2.0f * dy); //OK
	}
	for (int j = 1; j <= size_y - 2; j++) { //lado direito da rede principal ou seja i=(size_x-1)
		int right = (size_x - 1) + j * stride;
		GridPoint point(right, size_x, size_y, false);
		int westwest = right - 2;
		array_out_x[right] = (3.0f * (array_in[point.C]) - 4.0f * (array_in[point.W]) + 1.0f * (array_in[westwest])) /
		                     (2.0f * dx); //backwar difference
		array_out_y[right] = ((array_in[point.N]) - (array_in[point.S])) / (2.0f * dy);
	}

	int kp;
// i=0 j=0 forward x forward y
	kp = 0 + 0 * size_x;
	array_out_x[kp] = (-3.0f * (array_in[kp]) + 4.0f * (array_in[kp + 1]) - 1.0f * (array_in[kp + 2])) / (2.0f * dx);
	array_out_y[kp] = (-3.0f * (array_in[kp]) + 4.0f * (array_in[kp + 1 * size_x]) - 1.0f * (array_in[kp + 2 * size_x])) / (2.0f * dy);
// i=(size_x-1) j=0 backward x forward y
	kp = (size_x - 1) + 0 * size_x;
	array_out_x[kp] = (3.0f * (array_in[kp]) - 4.0f * (array_in[kp - 1]) + 1.0f * (array_in[kp - 2])) / (2.0f * dx);
	array_out_y[kp] = (-3.0f * (array_in[kp]) + 4.0f * (array_in[kp + size_x]) - 1.0f * (array_in[kp + 2 * size_x])) / (2.0f * dy);
// i=0 j=(size_y-1) forward x backward y
	kp = 0 + (size_y - 1) * size_x;
	array_out_x[kp] = (-3.0f * (array_in[kp]) + 4.0f * (array_in[kp + 1]) - 1.0f * (array_in[kp + 2])) / (2.0f * dx);
	array_out_y[kp] = (3.0f * (array_in[kp]) - 4.0f * (array_in[kp - size_x]) + 1.0f * (array_in[kp - 2 * size_x])) / (2.0f * dy);
// i=(size_x-1) j=(size_y-1) backward x backward y
	kp = (size_x - 1) + (size_y - 1) * size_x;
	array_out_x[kp] = (3.0f * (array_in[kp]) - 4.0f * (array_in[kp - 1]) + 1.0f * (array_in[kp - 2])) / (2.0f * dx);
	array_out_y[kp] = (3.0f * (array_in[kp]) - 4.0f * (array_in[kp - size_x]) + 1.0f * (array_in[kp - 2 * size_x])) / (2.0f * dy);
}


