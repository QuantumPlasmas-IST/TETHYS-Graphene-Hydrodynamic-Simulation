#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <iomanip>


#include <H5Cpp.h>

#include "TethysLib.h"
#include "Tethys1DLib.h"
#include "BoundaryLib.h"


#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace H5;
using namespace std;



#include <gsl/gsl_linalg.h>
const FloatType      HDF5FLOAT(PredType::NATIVE_FLOAT);
/*
 * Ficheiro de teste para parte difusiva 1D
 *
 * Implementação de Crank-Nicholson recorrendo à biblioteca GSL
 *
 * Depois testaremos em 2D
 *
*/

int main(int argc, char **argv){
	const int npoints=101; 							// number of spatial points
	float t=0.0;
	float dx=0.01f;								// spatial discretisation
	float dt=0.001f;								// time step

//	Fluid1D graph(npoints,40.0f,0.1f);

	GrapheneFluid1D	graph(npoints, 40.0f,10.0f, 0.1f, 0.0f);
	graph.CFLCondition();
	dx=graph.GetDx();
	dt=graph.GetDt();

	graph.InitialCondTest();



	int data_save_mode=1;

	float v_test[npoints];
	float mean=dx*npoints/2.0f;
	float sigma=dx*npoints/10.0f;
	for (int i = 0; i < npoints; i++ ){
		//v_test[i] = 0.4f*exp(-0.5f*((i*dx-mean)*(i*dx-mean)/(sigma*sigma)))/sigma;
		v_test[i] = 1.0f+tanh(20.0f*(dx*i-0.5));
		graph.Den[i]=1.0f;
	}


	float a = 0.1f*(0.5f*dt/(dx*dx));
	//versao pontos interiores
/*
    gsl_vector *x_solution = gsl_vector_alloc (npoints-2);
	gsl_vector *b_vector = gsl_vector_alloc (npoints-2);
	gsl_vector *diagonal = gsl_vector_alloc (npoints-2);
	gsl_vector *sub_diagonal = gsl_vector_alloc (npoints-2-1);
	gsl_vector *sup_diagonal = gsl_vector_alloc (npoints-2-1);
	gsl_vector_set_all(diagonal, 1.0f+2.0f*a);
	gsl_vector_set_all(sub_diagonal, -1.0f*a);
	gsl_vector_set_all(sup_diagonal, -1.0f*a);
*/
//versao com condicoes fronteira
	gsl_vector *x_solution = gsl_vector_alloc (npoints);
	gsl_vector *b_vector = gsl_vector_alloc (npoints);
	gsl_vector *diagonal = gsl_vector_alloc (npoints);
	gsl_vector *sub_diagonal = gsl_vector_alloc (npoints-1);
	gsl_vector *sup_diagonal = gsl_vector_alloc (npoints-1);
	gsl_vector_set_all(diagonal, 1.0f+2.0f*a);
	gsl_vector_set_all(sub_diagonal, -1.0f*a);
	gsl_vector_set_all(sup_diagonal, -1.0f*a);
	gsl_vector_set(diagonal, 0, 1.0f);
	gsl_vector_set(diagonal, npoints-1, 1.0f);
	gsl_vector_set(sup_diagonal, 0, 0.0f);
	gsl_vector_set(sub_diagonal, npoints-2, 0.0f);
	float aux=0.0f;


	int time_step=0;

	float t_max=2.0f;
	//Main cycle

	FILE * f = fopen ("testDiffusion1D.dat", "w");

	BoundaryCondition BC;

	while(t <= t_max ) // throw exception para nan / in
	{
		++time_step;
		t += dt;



//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
		for (int i = 1; i <= npoints - 2; i++) {
			aux = (1.0f - 2.0f * a) * v_test[i] + a * v_test[i - 1] + a * v_test[i + 1];
			gsl_vector_set(b_vector, i - 1, aux);
		}
*/
		for (int i = 1; i <= npoints - 2; i++) {
			aux = (1.0f - 2.0f * a) * graph.Vel[i] + a * graph.Vel[i-1] + a * graph.Vel[i+1];
			//aux = (1.0f - 2.0f * a) * v_test[i] + a * v_test[i - 1] + a * v_test[i + 1];
			gsl_vector_set(b_vector, i, aux);
		}
		//gsl_vector_set(b_vector, 0,  v_test[0]);
		//gsl_vector_set(b_vector, npoints-1,  v_test[npoints-1]);
		gsl_vector_set(b_vector, 0,   0.0f);
		gsl_vector_set(b_vector, npoints-1,   2.0f);
		//gsl_vector_set(b_vector, 0,   graph.Vel[0]);
		//gsl_vector_set(b_vector, npoints-1,   graph.Vel[npoints-1]);
		gsl_linalg_solve_tridiag(diagonal, sup_diagonal, sub_diagonal, b_vector, x_solution);
		gsl_vector_fprintf(f, x_solution, "%.4g");
		fprintf(f, "\n\n");
//		for (int i = 0; i <= npoints - 3; i++) {
//			v_test[i + 1] = gsl_vector_get(x_solution, i);
//		}

		for (int i = 0; i <= npoints - 1; i++) {
			graph.Vel[i] = gsl_vector_get(x_solution, i);
			//v_test[i] = gsl_vector_get(x_solution, i);
		}

		cout<<"um"<<graph.Vel[0]<<endl;
		graph.Richtmyer();
graph.Smooth(2);
		cout<<"dois"<<graph.Vel[0]<<endl;;
		BC.XFree(graph);
		cout<<"tres"<<graph.Vel[0]<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	gsl_vector_free (x_solution);
	gsl_vector_free (b_vector);
	gsl_vector_free (diagonal);
	gsl_vector_free (sub_diagonal);
	gsl_vector_free (sup_diagonal);
	fclose (f);
	return 0;
}




