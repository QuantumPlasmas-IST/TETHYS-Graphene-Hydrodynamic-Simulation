#ifndef TETHYSMATHLIB_H
#define TETHYSMATHLIB_H

#include <H5Cpp.h>
#include "TethysLib.h"

using namespace std;

void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);
constexpr float Gauss_Kernel(int position , float t); //
constexpr float Gauss_Kernel_Derivative(int position , float t); //
float Signal_Average(int n, float dt, const float * f);
float Integral_1_D(int n, float ds, const float * f);
float Integral_2_D(int n, int m, float dx, float dy, const float * f);
/* Average moving filter for the smoothing of 1D simulation, suppressing the spurious oscillations inherent to the 2nd order solver*/
void Average_Filter(const float * vec_in, float * vec_out, int size , int width );


#endif //TETHYSMATHLIB_H
