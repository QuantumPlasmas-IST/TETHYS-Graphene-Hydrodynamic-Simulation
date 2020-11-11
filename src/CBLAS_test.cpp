//
// Created by pcosme on 10/11/2020.
//

#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <string>
#include <random>
#include <exception>


//#include <gsl/gsl_cblas.h>

using namespace  std;

extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

//extern "C" {
//extern int sgemm_(char TRANSA, char TRANSB, int M, int N, int K, float  ALPHA,float* A , int LDA, float*B , int LDB, float BETA , float* C , int LDC );
//}
extern "C" {
extern int sgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, float*  ALPHA,float* A , int* LDA, float*B , int* LDB, float* BETA , float* C , int* LDC );
}



int main (){
	int lda = 2;
	float A[] = { 0.11, 0.12, 0.13,
	              0.21, 0.22, 0.23 };
	int ldb = 3;
	float B[] = { 1011, 1012,
	              1021, 1022,
	              1031, 1032 };
	int ldc = 2;
	float C[] = { 0.00, 0.00,
	              0.00, 0.00 };
	/* Compute C = A B */

	//cblas_sgemm (CblasRowMajor,CblasNoTrans, CblasNoTrans, 2, 2, 3,1.0, A, lda, B, ldb, 0.0, C, ldc);
	char Nchar='N';
	int M=2;
	int N=2;
	int K=3;
	float alpha=1.0;
	float beta=0.0;
	sgemm_(&Nchar, &Nchar, &M, &N, &K, &alpha, A ,&lda, B , &ldb, &beta , C , &ldc);

	printf ("[ %g, %g\n", C[0], C[1]);
	printf ("  %g, %g ]\n", C[2], C[3]);
/*
 * [ 367.76, 368.12
 *  674.06, 674.72 ]
  */


/*
	float * R;
	float * M;
	float * u;
	float * u_new;
	float dx=0.01, dt=0.01;

	int size=100;

	R = new float[size]();
	M = new float[size*size]();

	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++) {
			if(i==j){
				M[i+j*size] = -2;
			}
			if(i==j-1||i==j+1){
				M[i+j*size] = 1;
			}
		}
	}

	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++) {
			cout <<M[i+j*size]<<" ";
		}
		cout <<"\n";
	}


	u = new float[size]();
	u_new = new float[size]();

	for(int i=0;i<size;i++){
		float x=i*dx;
		if(x<=0.6&&x>=0.4) {
			R[i] = dx*dx*sqrt(-(x - .5) * (x - .5) + 0.01);
		}
	}


	ofstream myfile;
	myfile.open ("example_cblas.txt");

	for(float t=0;t<=1;t+=dt){

		//M u_new = R



		cblas_sgemv(CblasRowMajor,CblasNoTrans,size,size,1.0,M,100,R,1,0.0,u_new,1);
		cblas_scopy(size, u_new, 1, u, 1); //SCOPY copies a vector, x, to a vector, y.uses unrolled loops for increments equal to 1.

		for(int i=0;i<size;i++){
			myfile <<i*dx<<"\t";
			myfile <<u[i]<<"\n";
		}
		myfile << "\n\n";
	}
	myfile.close();

*/
	return 0;
}


