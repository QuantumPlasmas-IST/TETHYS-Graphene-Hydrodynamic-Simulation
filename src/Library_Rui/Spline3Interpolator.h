#ifndef __Spline3Interpolator__
#define __Spline3Interpolator__

#include "DataPoints.h"
#include "FCmatrixBanded.h"

class Spline3Interpolator : public DataPoints {
public:
	Spline3Interpolator(int fN = 0, double *fx = NULL, double *fy= NULL);
	Spline3Interpolator(std::vector<std::pair<double, double>> xy);
	Spline3Interpolator(const Spline3Interpolator&);

	~Spline3Interpolator();


	double Interpolate(double x);


private:
	void SetCurvatureLines();

	Vec K; //2nd derivatives
};

#endif