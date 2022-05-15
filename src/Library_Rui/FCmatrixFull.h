#ifndef __FCmatrixFull__
#define __FCmatrixFull__

#include "FCmatrix.h"
#include "FCmatrixSparse.h"
#include "FCmatrixAlgorithms.h"
#include "EqSolver.h"


class FCmatrixFull : public FCmatrix {
	
	public:

	//constructors
	FCmatrixFull(double ** fM, int fm, int fn);
	FCmatrixFull(double * fM, int fm, int fn);
	//default constructor
	FCmatrixFull(int fm = 3, int fn = 3, double a = 0);
	FCmatrixFull(std::vector<Vec>);

	//copy constructor
	FCmatrixFull(const FCmatrix&);

	//destructor
	~FCmatrixFull()=default;

	//operators
	void 		 operator=(const FCmatrixFull&);
	void 		 operator=(const FCmatrixSparse&);
	FCmatrixFull operator+(const FCmatrixFull&);
	FCmatrixFull operator-(const FCmatrixFull&);
	FCmatrixFull operator*(const FCmatrixFull&);
	FCmatrixFull operator*(double lambda);
	Vec 		 operator*(const Vec&);
	void         operator*=(double lambda);
	void 		 operator*=(const FCmatrixFull&);
	void 		 operator-=(const FCmatrixFull&);
	void 		 operator+=(const FCmatrixFull&);
	Vec& 		 operator[](int i);
	Vec 		 operator[](int i) const;

	//methods

	//virtual inherited
	Vec    GetRow(int i) 		const;
	Vec    GetCol(int i) 		const;
	int    GetRowN()     		const;
	int    GetColN()     		const;
	int    GetRowMax(int i = 0); //row max element index
	int    GetColMax(int j = 0);
	double Determinant();
	void   swapRows(int, int);
	void   swapCols(int, int);
	double   Norm1();
	double   NormInf();

	FCmatrixFull GetTranspose();
	FCmatrixFull GetInverse();

	std::vector<int> GetRowIndices() const;
	std::vector<int> GetColIndices() const;

	//print methods
	void Print() 	    const;
	void PrintIndices() const;

	//friend methods
	friend std::ostream& operator<<(std::ostream&, const FCmatrixFull&);

private:
	int *rowindices; //array of row indices (0,1,2,...)
	int *colindices;
};


#endif