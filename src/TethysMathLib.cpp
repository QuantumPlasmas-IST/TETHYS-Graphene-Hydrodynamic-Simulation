#include "TethysMathLib.h"

using namespace std;




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
	for(int k= 1 + n; k <= n * m - n - 2; k++) { //correr a grelha principal evitando as fronteiras
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
	avg = avg/(n * dt);
	return avg;
}
constexpr float Gauss_Kernel(int position , float t){
	return exp(-0.5f*position*position/t)/(sqrt(2.0f*MAT_PI*t));
}

constexpr float Gauss_Kernel_Derivative(int position , float t){
	return (-position*exp(-0.5f*position*position/t)/t)/(sqrt(2.0f*MAT_PI*t));
}

void Convolve_Gauss(int type, int m, float t, const float * in, float * out, int size){
	if(type==0){
		for(int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-m; k <= m; k++){
					out[i] += in[i-k] * Gauss_Kernel(k, t);
				}
			}
		}
	}
	if(type==1){
		for(int i=0;i<size;i++){
			if(i >= m && i < size - m){
			for(int k=-m; k <= m; k++){
					out[i] += in[i-k] * Gauss_Kernel_Derivative(k, t);
				}
			}
		}
	}
}
