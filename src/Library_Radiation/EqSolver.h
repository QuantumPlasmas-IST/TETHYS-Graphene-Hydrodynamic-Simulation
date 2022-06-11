#ifndef __EqSolver__
#define __EqSolver__

#include "FCmatrixFull.h"
#include "FCmatrixSparse.h"
#include "FCmatrixBanded.h"
#include "FCmatrixAlgorithms.h"

class EqSolver {
public:
	EqSolver();
	EqSolver(const FCmatrix&, const Vec&); // matriz M e vector de constantes B
	~EqSolver() = default;

	void SetConstants(const Vec&);
	void SetMatrix(const FCmatrix&);

	Vec GaussEliminationSolver();
	Vec LUdecompositionSolver();
	Vec GaussSeidelIterator(int  tolerance=1e-2);

private:
	FCmatrix *M; //matriz de coeffs
	Vec b; //vector de constantes
};

#endif