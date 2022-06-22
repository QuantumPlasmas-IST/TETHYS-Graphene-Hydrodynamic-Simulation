#ifndef __FCmatrixSparse__
#define __FCmatrixSparse__

#include "FCmatrix.h"
#include "Tools1.h"

class FCmatrixSparse : public FCmatrix {
public:

	//constructors
	FCmatrixSparse();
	FCmatrixSparse(double ** fM, int fm, int fn);
	FCmatrixSparse(double * fM, int fm, int fn);
	FCmatrixSparse(std::vector<Vec> v);
	FCmatrixSparse(const FCmatrixSparse&);

	//destructos
	~FCmatrixSparse() = default;

	//methods

	//virtual inherited
	Vec  	GetA() 		  		const; //vector with all the values of the !=0 entries
	Vec  	GetB() 		  		const; //vector with the indices of the first values != 0 in each row, relatively to A
	Vec  	GetC() 	      		const; //vector with all the column indices of the != values
	Vec  	GetRow(int i) 		const;
	Vec  	GetCol(int i) 		const;
	int  	GetRowN() 	  		const;
	int  	GetColN()     		const;
	int  	GetRowMax(int i);
	int     GetColMax(int j);
	double  Determinant();

	std::vector<int> GetRowIndices() const;
	std::vector<int> GetColIndices() const;

	void SetA(const Vec&);
	void SetB(const Vec&);
	void SetC(const Vec&);
	void SetN(const int);
	//operators

	void 		   operator=(const FCmatrixSparse&);
	FCmatrixSparse operator*(double lambda);
	Vec&		   operator[](int i);
	Vec 		   operator[](int i) const;
	FCmatrixSparse operator+(const FCmatrix&);
	FCmatrixSparse operator-(const FCmatrix&);
	void 		   operator-=(const FCmatrix&);
	void 		   operator+=(const FCmatrix&);
	/*FCmatrixSparse operator*(const FCmatrix&);
	Vec 		   operator*(const Vec&);
	void           operator*=(double lambda);
	void 		   operator*=(const FCmatrix&);
	*/

	//prints
	void  Print() const;
	//friend methods
	friend std::ostream& operator<<(std::ostream&, const FCmatrixSparse&);

private:
	int n;
};

#endif