#ifndef __Vec__
#define __Vec__

#include <iostream>
#include <cstdio>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <vector>

class Vec{
public:
	Vec(int i=0, double x=0);
	Vec(int i, const double* x);
	Vec(int i, const int* x);
	Vec(const Vec&); //copy constructor
	Vec(std::vector<double>);

	~Vec(); //destructor

	void SetEntries (int, double*);
	int * GetIndices();

	double  operator[](int i) const;
	double& operator[](int i);

	void    operator=(const Vec&); //A=B
	Vec     operator+(const Vec&); //C=A+B
	void    operator+=(const Vec&);
	Vec     operator-(const Vec&); //C=A-B
	Vec     operator-(); //C=-A
	void    operator-=(const Vec&);
	Vec     operator*(const Vec&) const; //C=A*B
	Vec     operator*(const double) const; //A = B*5
	void    operator*=(const Vec&); //B=A*A
	Vec     operator!(); //normalize vec

	double  dot(const Vec&); //double result = a.dot(b)
	Vec     ex(const Vec&); //C = a.ex(b) externo
	void    swap(int, int);
	int     size() const;
	double  mod() const; //norma 2
	double  AbsMax() const; //norma 1
	double  sumAbs(const Vec&);
	double* data();
	void    Print();

	friend std::ostream& operator<<(std::ostream&, const Vec&);
	friend Vec operator*(double, const Vec&); //5*A


private:
	int N;           //number of elements
	double* entries; //array
	int* indices; //indices
};

void swap(Vec&, Vec&);

#endif