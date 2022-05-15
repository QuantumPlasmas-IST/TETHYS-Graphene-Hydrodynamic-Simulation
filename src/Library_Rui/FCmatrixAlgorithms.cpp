#include "FCmatrixAlgorithms.h"

//#define DEBUG

void FCmatrixAlgorithms::LUdecomposition(FCmatrix& M_i, Vec& b, std::vector<int>& index){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	FCmatrixFull M(M_i);
	double lambda;
	int n = M[0].size();
	int j_aceite;

	//creating a matrix with nxn 0's
	FCmatrixFull L(n,n,0);

	//calculating scale factor for each line
	std::vector<double> scale_factor;
	for(int j=0; j<n; j++) //cicle on rows
		scale_factor.push_back(M[j][M.GetRowMax(j)]);

	for(int k=0; k<n-1; k++){ //line used as pivot
		j_aceite = k;
		for(int j=k; j<n; j++){  //testing other lines for pivot
			if(fabs(M[j][k]/scale_factor[j]) > fabs(M[j_aceite][k]/scale_factor[j_aceite])) //testing row with biggest relative size
				if(scale_factor[j_aceite]!= 0) j_aceite = j;
		}
		if(fabs(M[j_aceite][k]/scale_factor[j_aceite])<1e-3) //mal condicionado
			//throw std::invalid_argument(Form("[%s] Matrix not well conditioned!\n", __PRETTY_FUNCTION__));

		if(k != j_aceite){ //swaping the trade made above
			std::swap(scale_factor[k], scale_factor[j_aceite]);
			M.swapRows(k,j_aceite);
			L.swapRows(k, j_aceite);
			b.swap(k,j_aceite);
		}

		for(int i=k+1; i<n; i++){ //actual elimination phase
			lambda = M[i][k] / M[k][k];
			M[i] -= M[k] * lambda;
			L[i][k] = lambda;
		}
	}

	for(int i=0; i<n; i++)
		M_i[i] = M[i] + L[i];
}

int obi_1 = 0;

void FCmatrixAlgorithms::GaussElimination(FCmatrix& M_i, Vec& b, std::vector<int>& index){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	FCmatrixFull M(M_i);
	double lambda;
	int n = M[0].size();
	int j_aceite;
	obi_1 = 0;

	//calculating scale factor for each line
	std::vector<double> scale_factor;
	for(int j=0; j<n; j++) //cicle on rows
		scale_factor.push_back(M[j][M.GetRowMax(j)]);

	for(int k=0; k<n-1; k++){ //line used as pivot
		j_aceite = k;
		for(int j=k; j<n; j++){  //testing other lines for pivot
			if(fabs(M[j][k]/scale_factor[j]) > fabs(M[j_aceite][k]/scale_factor[j_aceite])) //testing row with biggest relative size
				if(scale_factor[j_aceite]!= 0) j_aceite = j;
		}
		if(fabs(M[j_aceite][k]/scale_factor[j_aceite])<1e-5) //mal condicionado
			obi_1++;
			//throw std::invalid_argument(Form("[%s] Matrix not well conditioned!\n", __PRETTY_FUNCTION__));
			//std::cout << "[Warning] Matrix is not well conditioned" << std::endl;

		if(k != j_aceite){ //swaping the trade made above
			std::swap(scale_factor[k], scale_factor[j_aceite]);
			M.swapRows(k,j_aceite);
			b.swap(k,j_aceite);
		}

		for(int i=k+1; i<n; i++){ //actual elimination phase
			lambda = M[i][k] / M[k][k];
			M[i] -= M[k] * lambda;
			b[i] -= b[k] * lambda;
		}
	}

	for(int i=0; i<n; i++){
		M_i[i] = M[i];
	}
}

int FCmatrixAlgorithms::GaussEliminationSimple(FCmatrix& M_i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	FCmatrixFull M(M_i);
	double lambda;
	int n = M[0].size();
	int j_aceite;
	int menos_mais = 1;

	//calculating scale factor for each line
	std::vector<double> scale_factor;
	for(int j=0; j<n; j++) //cicle on rows
		scale_factor.push_back(M[j][M.GetRowMax(j)]);

	for(int k=0; k<n-1; k++){ //line used as pivot
		j_aceite = k;
		for(int j=k; j<n; j++){  //testing other lines for pivot
			if(fabs(M[j][k]/scale_factor[j]) > fabs(M[j_aceite][k]/scale_factor[j_aceite])) //testing row with biggest relative size
				if(scale_factor[j_aceite]!= 0) j_aceite = j;
		}
		if(fabs(M[j_aceite][k]/scale_factor[j_aceite])<1e-3) //mal condicionado
			throw std::invalid_argument(Form("[%s] Matrix not well conditioned!\n", __PRETTY_FUNCTION__));

		if(k != j_aceite){ //swaping the trade made above
			std::swap(scale_factor[k], scale_factor[j_aceite]);
			M.swapRows(k,j_aceite);
			if(menos_mais==1) menos_mais = -1;
			else menos_mais = 1;
		}

		for(int i=k+1; i<n; i++){ //actual elimination phase
			lambda = M[i][k] / M[k][k];
			M[i] -= M[k] * lambda;
		}
	}

	for(int i=0; i<n; i++){
		M_i[i] = M[i];
	}

	return menos_mais;
}