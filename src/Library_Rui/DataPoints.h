#ifndef __DataPoints__
#define __DataPoints__

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "Tools1.h"

class DataPoints{
public:
	DataPoints();
	DataPoints(int, double*, double*);
	DataPoints(std::vector<std::pair<double, double>> xy);

	virtual ~DataPoints();

	virtual double Interpolate(double x);
	virtual void Print(std::string FILE="Files_Images_LET/datapoints.csv");

	double GetMinX(){ return xmin;};
	double GetMaxX(){ return xmax;};
	double GetMinY(){ return ymin;};
	double GetMaxY(){ return ymax;};

	std::vector<std::pair<double,double>> DerivativeVector(std::vector<std::pair<double,double>> v);
	std::vector<std::pair<double,double>> DerivativeVector(int M,std::vector<std::pair<double,double>> v);
	std::vector<std::pair<double,double>> GetMovingAverage(int M);

	double Integrate(double a=-1e27, double b=1e27); //trapezoidal para quaisqueres pontos

//protected:
	int N; // number of data points
	double* x; // arrays
	double* y; // arrays

	double xmin;
	double xmax;
	double ymin;
	double ymax;
};

#endif