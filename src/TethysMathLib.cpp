/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/TethysMathLib.h"

using namespace std;

float Sound_Velocity_Anisotropy(float x, float s) {
	return s;
}
float Sound_Velocity_Anisotropy(float x, float y, float s) {
	float s_mod;
	s_mod=s;
	return s_mod;
}

float Integral_1_D(int n, float ds, const float * f){
	float itg=0.0;
	for(int j=1; j < n / 2; j++){
		itg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	itg = itg*ds/3.0f;
	return itg;
}

float Integral_2_D(int n, int m, float dx, float dy, const float * f){
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


void Average_Filter(const float * vec_in, float * vec_out, int size , int width ){
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



float Signal_Average(int n, float dt, const float * f){
	float avg=0.0;
	for(int j=1; j < n / 2; j++){
		avg += f[2*j-2] + 4*f[2*j-1] + f[2*j];
	}
	avg = avg*dt/3.0f;
	avg = avg/(static_cast<float>(n) * dt);
	return avg;
}
float Gauss_Kernel(int position , float t){
	auto pos=static_cast<float>(position);
	return exp(-0.5f * pos * pos / t) / (sqrt(2.0f * MAT_PI * t));
}

float Gauss_Kernel_Derivative(int position , float t){
	auto pos=static_cast<float>(position);
	return (-pos * exp(-0.5f * pos * pos / t) / t) / (sqrt(2.0f * MAT_PI * t));
}

void Convolve_Gauss(unsigned int type, unsigned int m, float t, const float * in, float * out, unsigned long size){
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

float Stair_Case_Function(float x, float step_width, float smoothness) {
	return (0.5f*tanh(smoothness*((-(x - 1.0f)/step_width - floor(-(x - 1.0f)/step_width)) - 0.5f))/tanh(0.5f*smoothness) + 0.5f +floor(-(x - 1.0f)/step_width));
}
