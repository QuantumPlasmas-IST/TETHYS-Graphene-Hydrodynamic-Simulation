#ifndef __FCmatrixBanded__
#define __FCmatrixBanded__

#include "FCmatrix.h"

class FCmatrixBanded : public FCmatrix{

public:

  //constructors
  FCmatrixBanded();
  FCmatrixBanded(double **, int, int,int);
  FCmatrixBanded(double *, int, int, int);
  FCmatrixBanded(std::vector<Vec>&);
  FCmatrixBanded(const FCmatrixBanded&);
  FCmatrixBanded(const FCmatrix&);

  //destructor
  ~FCmatrixBanded() = default;

  //methods

  //virtual inherited
  double Determinant();
  int    GetRowN()              const;
  int    GetColN()              const;
  int    GetRowMax(int i=0);
  int    GetColMax(int j=0);
  int    GetUpperD()            const;
  int    GetLowerD()            const;
  Vec    GetRow(int)            const; 
  Vec    GetCol(int)            const;
  Vec    GetBand(int)           const;
  int    GetMainDiagonalIndex() const;

  std::vector<int> GetRowIndices() const;
  std::vector<int> GetColIndices() const;

  //operators
  Vec            operator[](int) const;
  Vec&           operator[](int);
  void           operator= (const FCmatrixBanded&);
  void           operator+=(const FCmatrixBanded&);
  void           operator-=(const FCmatrixBanded&);
  FCmatrixBanded operator+ (const FCmatrixBanded&);
  FCmatrixBanded operator- (const FCmatrixBanded&);
  FCmatrixBanded operator* (const FCmatrix&);
  FCmatrixBanded operator* (double lambda);

  //print methods
  void Print()      const;
  void PrintBands() const;

  //solver
  Vec Thomas(Vec) const; 

  //friend methods
  friend std::ostream& operator<<(std::ostream&, const FCmatrixBanded&);

private:
  int N;      //amt of numbers in the main diagonal
  int lowerD; //lowerD is the number of diagonals below the main diagonal
  int upperD; //upperD is the number of diagonals above the main diagonal
};

#endif
