#include "EqSolver.h"

//#define DEBUG

EqSolver::EqSolver() : b(), M(){;}

EqSolver::EqSolver(const FCmatrix& A, const Vec& vb) : b(vb){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	if(A.GetClassName().compare("FCmatrixFull") == 0)
		M = new FCmatrixFull(A);
	else if(A.GetClassName().compare("FCmatrixSparse") == 0)
		M = new FCmatrixFull((FCmatrixSparse&)A);
	else if(A.GetClassName().compare("FCmatrixBanded") == 0)
		M = new FCmatrixFull((FCmatrixBanded&)A);
	else
		throw std::invalid_argument("Matrix is not of any valid option for class\n");
}

void EqSolver::SetConstants(const Vec& vb){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	b = vb;
}

void EqSolver::SetMatrix(const FCmatrix& A){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	if(A.GetClassName().compare("FCmatrixFull") == 0){
		delete M;
		M = new FCmatrixFull((FCmatrixFull&)A);
	}
	if(A.GetClassName().compare("FCmatrixSparse") == 0){
		delete M;
		M = new FCmatrixFull((FCmatrixSparse&)A);
	}
	if(A.GetClassName().compare("FCmatrixBanded") == 0){
		delete M;
		M = new FCmatrixFull((FCmatrixBanded&)A);
	}
}

Vec EqSolver::GaussEliminationSolver(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	int n = M->GetRowN();
	if (n != M->GetColN()) throw std::invalid_argument("Matrix must be square for the usage of this method for specific solution!\n");
	Vec x((*M)[0].size());
	std::vector<int> index = M->GetRowIndices();

	FCmatrixAlgorithms::GaussElimination((*M), b, index);
	//M->Print();

	// backward solution (Ux=y)
	// loop on rows from end to begin
	for (int k=n-1; k>=0; k--){
		double sumX=0;

		// scan values of a row from diagonal to the right
		for (int j=k+1; j<n; j++)
			sumX += x[j] * (*M)[k][j];

		// compute solution values
		x[k] = (b[k] - sumX) / (*M)[k][k];
	}

	return x;
}

Vec EqSolver::LUdecompositionSolver(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	int n = M->GetRowN();
	if (n != M->GetColN()) throw std::invalid_argument("Matrix must be square for the usage of this method for specific solution!\n");
	Vec x((*M)[0].size());
	std::vector<int> index = M->GetRowIndices();

	FCmatrixAlgorithms::LUdecomposition((*M), b, index);


	// forward solution (Ly=b)
	//loop on rows
	std::vector<double> y(n);
	double sumC;
	for (int k=0; k<n; k++) {//loop on rows
		sumC = 0.;
		for (int i=0; i<k; i++)
			sumC += y[i] * (*M)[k][i];

		y[k] = b[k] - sumC;
	}

	// backward solution (Ux=y)
	double sumX=0;
	for (int k=n-1; k>=0; k--){
		sumX=0;

		// scan values of a row from diagonal to the right
		for (int j=k+1; j<n; j++)
			sumX += x[j] * (*M)[k][j];

		// compute solution values
		x[k] = (y[k] - sumX) / (*M)[k][k];
	}

	return x;
}

Vec EqSolver::GaussSeidelIterator(int tolerance){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	int n = M->GetRowN();
	if (n != M->GetColN()) 
		throw std::invalid_argument("Matrix must be square for the usage of this method for specific solution!\n");

	bool diagonal_dominant = 1;
	double sum_ndiag;
	for(int i=0; i<n; i++){
		sum_ndiag = 0;
		for(int j=0; j<n; j++){
			if(i!=j) 
				sum_ndiag += fabs((*M)[i][j]);
		}
		if(sum_ndiag >= fabs((*M)[i][i]) && diagonal_dominant ==1) 
			diagonal_dominant = 0;
	}

	if(diagonal_dominant == 0){
		for(int i=0; i<n; i++){
					sum_ndiag = 0;
			for(int j=0; j<n; j++){
				if(i!=j) 
					sum_ndiag += fabs((*M)[j][i]);
			}
			if(sum_ndiag >= fabs((*M)[i][i]) && diagonal_dominant ==1) 
				diagonal_dominant = 0;
		}
	}

	if(diagonal_dominant == 0)
		std::cout << "Matrix is not diagonal dominant neither by row nor by colums, so there is an high chance of this method not converging" << std::endl;

	Vec x(n);
	Vec x_km1(n);

	Vec x_aux(n); //zero's
	bool btol = false;
	int it = 0.;

	while (!btol && (it++ < 10)){
		x_aux = x;
		for (int i=0; i<n; i++) {
			x[i] = 0.;
			for (int j=0; j<n; j++){
				if (i != j) x[i] -= (*M)[i][j]*x[j];
			}

			x[i] += b[i];
			x[i] /= (*M)[i][i];

			//guarantee that all vector entries are converging equally
		}
		if ((x-x_aux).mod() < tolerance) btol = true;
		else btol = false;
		if(it==8) x_km1 = x;
	}

	double delta_x = (x_km1-x).mod();

	while (!btol && (it++ < 20)){
		x_aux = x;
		for (int i=0; i<n; i++) {
			x[i] = 0.;
			for (int j=0; j<n; j++)
				if (i != j) x[i] -= (*M)[i][j]*x[j];

			x[i] += b[i];
			x[i] /= (*M)[i][i];
		}
		//guarantee that vector is converging
		if ((x-x_aux).mod() < tolerance) btol = true;
		else btol = false;
		if(it==18) x_km1 = x;
	}
	double delta_xp = (x_km1-x).mod();

	//std::cout << delta_x/delta_xp << "x" << delta_x << "xp" << delta_xp << std::endl;

	double omega = 2 / (1 + sqrt(1-pow(delta_xp/delta_x,0.1)));
	//std::cout << "Omega_relaxation = "<< omega << std::endl;

	double x_atual;
	if(omega < 1.85){
		while (!btol && (it++ < 500)){
			x_aux = x;
			for (int i=0; i<n; i++) {
				x_atual = 0.;
				for (int j=0; j<n; j++)
					if (i != j) x_atual -= (*M)[i][j]*x[j];

				x_atual += b[i];
				x_atual /= (*M)[i][i];
				x[i] = omega*x_atual + (1-omega)*x[i];
			}
			//guarantee that all vector entries are converging equally
			if ((x-x_aux).mod() < tolerance) btol = true;
			else btol = false;
		}
	}

	return x;
}