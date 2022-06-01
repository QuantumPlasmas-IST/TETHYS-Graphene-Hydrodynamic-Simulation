#include "DataPoints.h"

//#define DEBUG

DataPoints::DataPoints() : N(1){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	x = new double[N];
	y = new double[N];
}

DataPoints::DataPoints(int N_i, double* x_in, double* y_in) : N(N_i) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(N<0)
		throw std::invalid_argument("You can't have less than 0 data points\n");

	x = new double[N];
	y = new double[N];

	for(int i=0; i<N; i++){
		x[i] = x_in[i];
		y[i] = y_in[i];
	}

	xmin = x[0];
	xmax = x[0];

	for(int i=0; i<N; i++){
		if(x[i] < xmin) xmin = x[i];
		if(x[i] > xmax) xmax = x[i];
	}

	ymin = x[0];
	ymax = x[0];

	for(int i=0; i<N; i++){
		if(y[i] < ymin) ymin = y[i];
		if(y[i] > ymax) ymax = y[i];
	}

	sort2(N,x,y);
}

DataPoints::DataPoints(std::vector<std::pair<double, double>> xy) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	N = xy.size();

	if(N<0)
		throw std::invalid_argument("You can't have less than 0 data points\n");

	x = new double[N];
	y = new double[N];

	for(int i=0; i<N; i++){
		x[i] = xy[i].first;
		y[i] = xy[i].second;
	}

	xmin = x[0];
	xmax = x[0];

	for(int i=0; i<N; i++){
		if(x[i] < xmin) xmin = x[i];
		if(x[i] > xmax) xmax = x[i];
	}

	ymin = x[0];
	ymax = x[0];

	for(int i=0; i<N; i++){
		if(y[i] < ymin) ymin = y[i];
		if(y[i] > ymax) ymax = y[i];
	}

	sort2(N,x,y);
}

DataPoints::~DataPoints(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	delete x;
	delete y;
}

double DataPoints::Interpolate(double x){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	std::cout << "No interpolation defined with no method (define object in an interpolation class)" << std::endl;
	return 0.;
}

void DataPoints::Print(std::string FILE){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	const char* F = FILE.c_str(); 
	std::ofstream fil(F, std::fstream::out);
	for(int i=0; i<N; i++)
		fil << x[i] << "," << y[i] << std::endl;

	fil.close();
}

std::vector<std::pair<double,double>> DataPoints::GetMovingAverage(int M){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    std::vector<std::pair<double,double>> vec_f;
    double sum = 0.;
    double sumx = 0;

    for(int j =  M/2; j < N-M; j++) {
        sum = 0.;
        sumx = 0.;
        for (int i = j; i < j+M; ++i){
        	sum += y[i];
        	sumx += x[i];
        }
        if(j == M/2) for (int i = 0; i < M/2; i++) vec_f.push_back(std::make_pair(x[i], sum/M));
        vec_f.push_back(std::make_pair(sumx/M, sum/M));
        if(j == N-M-1) for (int i = N-M; i < N; i++) vec_f.push_back(std::make_pair(x[i], sum/M));
    }
    
    return vec_f;
}

std::vector<std::pair<double,double>> DataPoints::DerivativeVector(int M,std::vector<std::pair<double,double>> v){
    std::pair<double,double> temp_pair;
    std::vector<std::pair<double,double>> firstDer;
    double temp = 0;

    for(int j =  M/2+1; j < N-M/2-1; j++){
    	temp = (v[j+1].second - v[j-1].second)/(v[j+1].first-v[j-1].first);
        temp_pair = std::make_pair(v[j].first,temp);
        if(j == M/2+1) for (int i = 0; i < M/2+1; i++) firstDer.push_back(temp_pair);
        firstDer.push_back(temp_pair);
        if(j == N-M/2-2) for (int i = N-M/2-1; i < N; i++) firstDer.push_back(temp_pair);
    }

    return firstDer;
}

std::vector<std::pair<double,double>> DataPoints::DerivativeVector(std::vector<std::pair<double,double>> v){
    std::pair<double,double> temp_pair;
    std::vector<std::pair<double,double>> firstDer;

    for(int i=1; i<N-1; i++){
        temp_pair = std::make_pair(v[i].first,(v[i+1].second - v[i-1].second)/(v[i+1].first-v[i-1].first));
        firstDer.push_back(temp_pair);
        std::cout << (v[i+1].second - v[i-1].second)/(v[i+1].first-v[i-1].first) << std::endl;
        if(i==1){
        	temp_pair = std::make_pair(v[0].first,(v[i+1].second - v[i-1].second)/(v[i+1].first-v[i-1].first));
        	firstDer.push_back(temp_pair);
        }
        if(i==N-2){
        	temp_pair = std::make_pair(v[N-1].first,(v[i+1].second - v[i-1].second)/(v[i+1].first-v[i-1].first));
        	firstDer.push_back(temp_pair);
        }
    }

    return firstDer;
}

double DataPoints::Integrate(double a, double b){
	double integral = 0;

	for(int i=0; i<N-1; i++){
		if((x[i]>a) && (x[i]<b)) integral += (y[i+1]+y[i])*(x[i+1]-x[i])/2;
		//cout <<  << << endl;
	}

	return integral;
}