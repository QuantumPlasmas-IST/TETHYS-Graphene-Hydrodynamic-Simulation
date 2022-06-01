#ifndef __FCmatrix__
#define __FCmatrix__

#include "Vec.h"
#include <vector>
#include <iostream>
#include <stdexcept>

class FCmatrix{
public:

	FCmatrix(); // flag = 0
	FCmatrix(double **, int, int); //rows, columns flag = 1
	FCmatrix(double *, int, int); //rows, columns flag = 2
	FCmatrix(const std::vector<Vec>&); // flag = 3
	FCmatrix(const FCmatrix&);
	~FCmatrix() = default;

	//method
	int GetSize() const;
	std::string GetClassName() const;
	virtual int GetRowN() const = 0; //# rows
	virtual int GetColN() const = 0; //# cols
	virtual Vec GetRow(int i) const= 0;
	virtual Vec GetCol(int i) const = 0;
	virtual double Determinant() = 0;
	virtual void Print() const; //print matrix M
	virtual int GetRowMax(int i = 0) = 0; //row max element index
	virtual int GetColMax(int j = 0) = 0;
	virtual std::vector<int> GetRowIndices() const = 0;
	virtual std::vector<int> GetColIndices() const = 0;

	//operators
	virtual Vec& operator[](int i) = 0;
	virtual Vec operator[](int i) const = 0;

protected:
	std::string classname;
	std::vector<Vec> M3;
};

#endif