#include <cstdlib>
#include <iostream>
#include "LaxWen2D.h"
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;

void WaterFluxX(float* input, float* output)
{
	output[0] = input[1];
	output[1] = input[1]*input[1]/input[0];
	output[2] = input[1]*input[2]/input[0];
}

void WaterFluxY(float* input, float* output)
{
	output[0] = input[2];
	output[1] = input[2]*input[1]/input[0];
	output[2] = input[2]*input[2]/input[0];
}

int main()
{
	// Basic Testing
	LaxWen2D Solver;
	Solver.DumpVariables();

	Solver.SetGrid(0.05, 501, 0.001, 101);
	Solver.SetTimeStep(0.000001);
	float* vec1 = new float[101];
	for(int i=0; i<101; i++)
		vec1[i] = i;
	Solver.SetBV(27, vec1);
	Solver.SetBV(2, vec1);
	Solver.DumpVariables();
	//Solver.DumpBV(2);
	//Solver.DumpBV(32);

	bool* vec_bool = new bool[12];
	for(int i=0; i<12; i++)
	{
		if(i%2 == 0)
			vec_bool[i] = true;
		else
			vec_bool[i] = false;
	}
	Solver.SetBV(vec_bool);
	Solver.DumpVariables();
	//Solver.DumpBV(2);
	//Solver.DumpBV(33);

	// Physical testing
	for(int i=0; i<12; i++)
		vec_bool[i] = true;
	LaxWen2D PhyTest(0.01, 201, 0.01, 201, 0.00001, vec_bool);
	PhyTest.DumpVariables();

	srand(time(NULL));
	float*** var_cur;
	var_cur = new float**[201];
	for (int i = 0; i < 201; i++)
	{
		var_cur[i] = new float*[201];
		for (int j = 0; j < 201; j++)
		{
			var_cur[i][j] = new float[3];
			for(int k=0; k<2; k++)
				var_cur[i][j][k] = 1;
			var_cur[i][j][2] = -0.005;
			for(int k=0; k<3; k++)
				var_cur[i][j][k] += 0.01*(float)rand()/(float)RAND_MAX;
		}
	}

	ofstream file_x0;
	file_x0.open("data/random_x0.txt");
	ofstream file_xM;
	file_xM.open("data/random_xM.txt");
	ofstream file_xL;
	file_xL.open("data/random_xL.txt");
	ofstream file_y0;
	file_y0.open("data/random_y0.txt");
	ofstream file_yM;
	file_yM.open("data/random_yM.txt");
	ofstream file_yL;
	file_yL.open("data/random_yL.txt");

	for(int m=0; m<100; m++)
	{
		for(int j=0; j<201; j++)
		{
			file_x0
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[0][j][0] 
			<< "\t" << var_cur[0][j][1]
			<< "\t" << var_cur[0][j][2] 
			<< endl;
			file_xM
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[100][j][0] 
			<< "\t" << var_cur[100][j][1]
			<< "\t" << var_cur[100][j][2] 
			<< endl;
			file_xL
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[200][j][0] 
			<< "\t" << var_cur[200][j][1]
			<< "\t" << var_cur[200][j][2] 
			<< endl;
		}
		for(int i=0; i<201; i++)
		{
			file_y0
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][0][0] 
			<< "\t" << var_cur[i][0][1]
			<< "\t" << var_cur[i][0][2] 
			<< endl;
			file_yM
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][100][0] 
			<< "\t" << var_cur[i][100][1]
			<< "\t" << var_cur[i][100][2] 
			<< endl;
			file_yL
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][200][0] 
			<< "\t" << var_cur[i][200][1]
			<< "\t" << var_cur[i][200][2] 
			<< endl;
		}
		PhyTest.Rich2(var_cur, WaterFluxX, WaterFluxY);
	}
	file_x0.close();
	file_xM.close();
	file_xL.close();
	file_y0.close();
	file_yM.close();
	file_yL.close();

	// Boundary conditions different 
	file_x0.open("data/sin_x0.txt");
	file_xM.open("data/sin_xM.txt");
	file_xL.open("data/sin_xL.txt");
	file_y0.open("data/sin_y0.txt");
	file_yM.open("data/sin_yM.txt");
	file_yL.open("data/sin_yL.txt");

	vec_bool[0] = false;
	vec_bool[1] = false;
	vec_bool[8] = false;
	vec_bool[11] = false;
	PhyTest.SetBV(vec_bool);

	float* vec_bv = new float[201];
	for (int i=0; i < 201; i++)
		vec_bv[i] = 0;
	PhyTest.SetBV(22, vec_bv);
	PhyTest.SetBV(32, vec_bv);
	for (int i=0; i < 201; i++)
		vec_bv[i] = 1;
	
	PhyTest.DumpVariables();

	for(int m=0; m<10000; m++)
	{
		for(int j=0; j<201; j++)
			vec_bv[j] = 1+ 0.5*sin(0.001*m);
		PhyTest.SetBV(0, vec_bv);
		PhyTest.SetBV(1, vec_bv);

		for(int j=0; j<201; j++)
		{
			file_x0
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[0][j][0] 
			<< "\t" << var_cur[0][j][1]
			<< "\t" << var_cur[0][j][2] 
			<< endl;
			file_xM
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[100][j][0] 
			<< "\t" << var_cur[100][j][1]
			<< "\t" << var_cur[100][j][2] 
			<< endl;
			file_xL
			<< m*0.0001 << "\t" 
			<< j*0.01 << "\t" 
			<< var_cur[200][j][0] 
			<< "\t" << var_cur[200][j][1]
			<< "\t" << var_cur[200][j][2] 
			<< endl;
		}
		for(int i=0; i<201; i++)
		{
			file_y0
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][0][0] 
			<< "\t" << var_cur[i][0][1]
			<< "\t" << var_cur[i][0][2] 
			<< endl;
			file_yM
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][100][0] 
			<< "\t" << var_cur[i][100][1]
			<< "\t" << var_cur[i][100][2] 
			<< endl;
			file_yL
			<< m*0.0001 << "\t" 
			<< i*0.01 << "\t" 
			<< var_cur[i][200][0] 
			<< "\t" << var_cur[i][200][1]
			<< "\t" << var_cur[i][200][2] 
			<< endl;
		}
		file_x0 << endl;
		file_xM << endl;
		file_xL << endl;
		file_y0 << endl;
		file_yM << endl;
		file_yL << endl;
		
		PhyTest.Rich2(var_cur, WaterFluxX, WaterFluxY);
		if(m%1000 == 0)
			cout << "Progress Check... Now at " << m/100 << "% var_cur[100][100][1] = "<< var_cur[100][100][1] << endl; 

	}
	file_x0.close();
	file_xM.close();
	file_xL.close();
	file_y0.close();
	file_yM.close();
	file_yL.close();

	delete[]vec1;
	delete[]vec_bool;
	for(int i=0; i<201; i++)
	{
		for(int j=0; j<201; j++)
			delete[] var_cur[i][j];
		delete[] var_cur[i];
	}
	delete[]var_cur;
	delete[]vec_bv;

	return 0;
}