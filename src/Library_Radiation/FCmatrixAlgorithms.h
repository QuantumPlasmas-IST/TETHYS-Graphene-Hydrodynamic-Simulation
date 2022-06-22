#ifndef __FCmatrixAlgorithms__
#define __FCmatrixAlgorithms__

#include "FCmatrixFull.h"
#include "FCmatrixSparse.h"
#include "FCmatrixBanded.h"


class FCmatrixAlgorithms{
public:
    //decomposição LU com |L|=1
    static void LUdecomposition(FCmatrix&, Vec&, std::vector<int>&);

    //return triangular matrix and changed vector of constants
    static void GaussElimination(FCmatrix&, Vec&, std::vector<int>&);

    //return triangular matrix
    static int GaussEliminationSimple(FCmatrix& M_i);
};

extern int obi_1;

#endif