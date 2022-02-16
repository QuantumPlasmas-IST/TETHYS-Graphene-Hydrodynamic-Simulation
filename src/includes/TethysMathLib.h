/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

/*!@file
 * @brief Header file for general mathematical functions
 */

#ifndef TETHYSMATHLIB_H
#define TETHYSMATHLIB_H



#include <cmath>
#include <functional>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <H5Cpp.h>

using namespace std;


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979f
#endif

#ifndef MAT_EULER
#	define MAT_EULER 2.71828182845905f
#endif

class MathUtils{

protected:
	gsl_matrix * FDmatrix2 ; //matrices to perform finite differenceing
	gsl_matrix * FDmatrix3 ;
	void SetFDmatrix2(int size);
	void SetFDmatrix3(int size);

public:
	MathUtils()=default;
	~MathUtils()=default;


/*!
 * @brief 1D Function of local anisotropy of S
 *
 * Function to implement the spatial variation of the sound velocity S(x) in 1D
 * corresponding to a variation of substrate permittivity or even the description of a multi gated system.
 *
 * @param x position along x
 * @param s nominal sound velocity
 *
 * @return S(x) the value of local sound velocity
 * */
static float Sound_Velocity_Anisotropy(float x, float s);

/*!
 * @brief 2D Function of local anisotropy of S
 *
 * Function to implement the spatial variation of the sound velocity S(x,y) in 2D
 * corresponding to a variation of substrate permittivity or even the description of a multi gated system.
 *
 * @param x position along x
 * @param y position along y
 * @param s nominal sound velocity
 *
 * @return S(x,y) the value of local sound velocity

 * */
static float Sound_Velocity_Anisotropy(float x, float y, float s);

/*!
 * @brief Smooth staircase function
 *
 *
 * This function is defined as @f[ f(x) = \frac{\tanh a(-\frac{x-1}{w} - \lfloor-\frac{x-1}{w})\rfloor - 1/2)}{2\tanh a/2 } + 1/2 +\lfloor-\frac{x-1}{w}\rfloor @f]
 * corresponding to a smooth descending staircase.
 *
 * @param x position
 * @param step_width width of each step, @f$ w @f$
 * @param smoothness controls the smoothness of the jump,  @f$ a @f$
 *
 * @return @f$ f(x) @f$
 *
 * */
static float Stair_Case_Function(float x, float step_width, float smoothness);

/*!
 * @brief Performs 1D convolution of array with a gaussian kernel or its derivative
 *
 *
 * @param type if type=0 convolve gaussian kernel, if type=1 convolve the first derivative of gaussian kernel
 * @param m window filter width
 * @param t standard deviation of gaussian kernel
 * @param in input array
 * @param out output array where to save the convolution
 * @param size array size
 *
 * @see Gauss_Kernel
 * @see Gauss_Kernel_Derivative
 * */
static void Convolve_Gauss(unsigned int type, unsigned int m, float t, const float * in, float * out, unsigned long size);
//void Convolve_Gauss(int type, float m, float t, float * in, float * out, int size);

/*!
 * @brief Gaussian kernel function
 *
 * This function implements a gaussian kernel defined as @f[ g(x) = \frac{1}{\sqrt{2\pi t}}e^{-\frac{x^2}{2t}} @f]
 *
 * @param pos value of variable @f$ x @f$
 * @param t standard deviation of gaussian kernel
 *
 * @return @f$ g(x) @f$
 * */
static float Gauss_Kernel(int pos , float t); //

/*!
 * @brief Gaussian kernel function derivative
 *
 * This function implements the derivative of a gaussian kernel defined as @f[ g'(x) = -\frac{x}{t\sqrt{2\pi t}}e^{-\frac{x^2}{2t}} @f]
 *
 * @param pos value of variable @f$ x @f$
 * @param t standard deviation of gaussian kernel
 *
 * @return @f$ g'(x) @f$
 * */
static float Gauss_Kernel_Derivative(int pos , float t); //


/*!
 * @brief Time average of time series or signal
 *
 * Used to perform time average on a set om sampled values of a signal or function
 *
 * @param n Number of points to integrate
 * @param dt time between points
 * @param f Pointer to array to average
 *
 * @return Time average of sampled function
 * */
static float Signal_Average(int n, float dt, const float * f);

/*!
 * @brief 1D Numerical integration using composite Simpson's rule
 *
 * @param n Number of points to integrate
 * @param ds Step between points
 * @param f Pointer to array of integrand
 *
 * @returns 4th order approximation of integral
 * */
static float Integral_1_D(int n, float ds, const float * f);


/*!
 * @brief 2D Numerical integration using composite Simpson's rule
 *
 * @param n Number of points along x to integrate
 * @param m Number of points along y to integrate
 * @param dx Step between points along x
 * @param dy Step between points along y
 * @param f Pointer to array of integrand
 *
 * @returns 4th order approximation of integral
 * */
static float Integral_2_D(int n, int m, float dx, float dy, const float * f);


/*!
 * @brief Moving average filter
 *
 * Average moving filter for the smoothing of 1D simulation, suppressing the spurious oscillations inherent to the 2nd order solver
 *
 * @param array_in Input array to be smoothed
 * @param array_out Output array to store the smoothed result
 * @param size Size of the arrays
 * @param width Window filter width
 *
 * */
static void Average_Filter(const float * array_in, float * array_out, int size , int width );


static float Signum(float x);

static void LaplacianField(const float *array_in, float *array_out, float dx, int size_x, int size_y);

static void GradientField(const float *array_in, float *array_out_x, float *array_out_y, float dx, float dy, int size_x, int size_y);
static void GradientField(const float *array_in, float *array_out, float dx, int size_x);
};

#endif