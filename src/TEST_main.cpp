#include "includes/SetUpParametersLib.h"
#include "includes/Fluid1DLib.h"
#include "includes/GrapheneFluid1DLib.h"

#include <functional>

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif

using namespace std;


float NonlinearWaveDESolver(float x, float dx, const float* Ic, const float* Par); // 4th order Runge-Kutta nonlinear wave DESolver


int main(int argc, char **argv) {
	float t = 0.0f;
	float dt;		// time step


	SetUpParameters parameters(20.0f,10.0f,0.0f,0.0f,0.0f,0.0f,0.0f,1,1.0);
	GrapheneFluid1D graph(parameters);


	graph.CflCondition();
	dt = graph.GetDt();

	graph.CreateFluidFile();
	graph.CreateHdf5File();


	/*
	// simplified B = 0.005f
	float offset = 1.000f;
	float A      = 0.010f;
	float f      = 4.000f;
	float phi    = MAT_PI;
	float c      = 21.1285f;
	*/

	/*
	// simplified B = 0.005f
	float offset = 1.0005f;
	float A      = 0.0106f;
	float f      = 2.0000f;
	float phi    = MAT_PI;
	float c      = 21.1943f;
	*/

	// simplified B = 0.005f
	float offset =  1.0000f;
	float A      =  0.3700f;
	float mu     =  0.5000f;
	float sigma  =  0.0550f;
	float c      = -21.2318f;

	// nonlinear wave approximate IC
	//std::function<float(float)> fden = [=](float x) { return offset + A*cos(2*MAT_PI*f*x+phi); };

	// soliton wave approximate IC
	std::function<float(float)> fden = [=](float x) { return offset + A/cosh((x-mu)/sigma); };
	std::function<float(float)> fvx  = [=](float x) { return c*(fden(x) - 1); };

	graph.InitialCondGeneral(fden, fvx);


	std::cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;


	graph.SetTmax(1.0f);
	float Tmax = graph.GetTmax();

	graph.CopyFields();
	graph.SaveSnapShot();
	while(t < Tmax) {

		t += dt;
		GrapheneFluid1D::TimeStepCounter++;

		graph.RungeKuttaTVD();

		if (parameters.SaveMode && graph.Snapshot()) {

			graph.CopyFields();
			graph.SaveSnapShot();
		}
		graph.WriteFluidFile(t);
	}
	if(parameters.SaveMode) {

		graph.WriteAttributes();
	}
	graph.CloseHdf5File();


	std::cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m" << endl;


	return 0;
}