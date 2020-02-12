#include<iostream>

using namespace std;

void TestFunc(float* input, float* output)
{
	for(int i = 0; i<5; i++)
		output[i] = 2* input[i];
}

int main()
{
	float* vec1 = new float[4];
	for(int i =0; i<4; i++)
		vec1[i] = i;
	cout << endl << "vec 1 components \n\n";
	for(int i =0; i<4; i++)
		cout << vec1[i] << endl;

	float* vec2 = new float[5];
	TestFunc(vec1, vec2);
	cout << endl << "vec 2 components \n\n";
	for(int i =0; i<4; i++)
		cout << vec2[i] << endl;

	float* vec3 = new float[4];
	TestFunc(vec2, vec3);
	cout << endl << "vec 2 components \n\n";
	for(int i =0; i<5; i++)
		cout << vec2[i] << endl;
	cout << endl << "vec 3 components \n\n";
	for(int i =0; i<4; i++)
		cout << vec3[i] << endl;

	float** mat;
	mat = new float* [3];

	for(int i=0; i<3; i++)
	{
		mat[i] = new float[5];
		for (int j = 0; j < 5; j++)
		{
			mat[i][j] = i*j;
		}
	}
	TestFunc(mat[0], mat[2]);

	for(int i = 0; i<3; i++)
	{
		cout << endl;
		for (int j = 0; j < 5; j++)
			cout << mat[i][j] << " ";

	}

	cout << endl;

	delete[] vec1;
	delete[] vec2;
	delete[] vec3;
	for(int i=0; i<3; i++)
		delete[] mat[i];
	delete[] mat;

	return 0;

}