#include "Spline3Interpolator.h"
#include <time.h>
#include "EqSolver.h"

//#define DEBUG


Spline3Interpolator::Spline3Interpolator(int fN, double *fx, double *fy) : DataPoints(fN,fx,fy), K(N) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	sort2(N,x,y);

	SetCurvatureLines(); //define segment interpolators
}

Spline3Interpolator::Spline3Interpolator(std::vector<std::pair<double, double>> xy) : DataPoints(xy), K(N) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif


	sort2(N,x,y);

	SetCurvatureLines(); //define segment interpolators
}

Spline3Interpolator::Spline3Interpolator(const Spline3Interpolator& S3) : Spline3Interpolator(S3.N, S3.x, S3.y){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
}


Spline3Interpolator::~Spline3Interpolator(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
}



void Spline3Interpolator::SetCurvatureLines() {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(N == 3){
		K[0] = 0;
		K[1] = 6*((y[0] - y[1]) / (x[0] - x[1]) - (y[1] - y[2]) / (x[1] - x[2])) / (2 * (x[0] - x[2]));
		K[2] = 0;
	}else{
		std::vector<Vec> bands =  {Vec(N-3), Vec(N-2), Vec(N-3)};
		
		bands[1][0] = 2*(x[0] - x[2]);

		for(int i=1; i< N-2; i++){
			bands[0][i-1] = x[i] - x[i+1];
			bands[1][i] = 2*(x[i] - x[i+2]);
		}
		bands[2] = bands[0];

		FCmatrixBanded MatK(bands);
		FCmatrixFull MF(MatK);

		Vec b(N-2);
		for(int i=0; i<N-2; i++){
			b[i] = 6*((y[i] - y[i+1]) / (x[i] - x[i+1]) -
					    (y[i+1] - y[i+2]) / (x[i+1] - x[i+2]));
		}

		Vec dumb(N-2);
		dumb = MatK.Thomas(b);

		for(int i=1; i<N-1; i++)
			K[i] = dumb[i-1];
	}
}

double Spline3Interpolator::Interpolate(double fx) {
	#ifdef DEBUG
	  printf("[%s] [x=%f]\n", __PRETTY_FUNCTION__,fx );
	#endif
	double declive1;
	double declive2;
	double b1;
	double b2;

	if(fx > xmax || fx < xmin){
		//std::cout << "(We advise against extrapolation)" << std::endl;
		declive1 = (y[1]-y[0])/(x[1]-x[0]);
		declive2 = (y[N-1]-y[N-2])/(x[N-1]-x[N-2]);
		b1 = y[0] - declive1*x[0];
		b2 = y[N-1] - declive2*x[N-1];
	}
	if(fx > xmax)
		return declive2*fx + b2;
	else if (fx < xmin)
		return declive1*fx + b1;

  	int i = 0;

  	for(int j=0; j<N-1; j++){
		if(fx > x[j] && fx <= x[j+1]){
      		i = j;
      		break;
    	}
  	}

  	//used for extrapolation to the right (0 used on the left)
  	if(fx>x[N-1]) i = N-2;

  	return  (K[i]/6)*((pow(fx-x[i+1],3))/(x[i]-x[i+1])-(fx-x[i+1])*(x[i]-x[i+1]))-
			(K[i+1]/6)*((pow(fx-x[i], 3))/(x[i]-x[i+1])-(fx-x[i])*(x[i]-x[i+1]))+
			(y[i]*(fx-x[i+1])-y[i+1]*(fx-x[i]))/(x[i]-x[i+1]);
}