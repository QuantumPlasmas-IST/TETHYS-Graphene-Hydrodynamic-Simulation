#include "FCmatrixFull.h"

//#define DEBUG


//-----------------------------------------------
//-----------------CONSTRUCTORS------------------
//-----------------------------------------------

FCmatrixFull::FCmatrixFull(double ** fM, int fm, int fn) : FCmatrix(fM, fm, fn){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!fM)
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));
	rowindices = new int[fm];
	colindices = new int[fn];
	for (int i = 0; i < fm; ++i)
		rowindices[i] = i;
	for (int i = 0; i < fn; ++i)
		colindices[i] = i;
	classname = "FCmatrixFull";
}

FCmatrixFull::FCmatrixFull(double * fM, int fm, int fn) : FCmatrix(fM, fm, fn) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!fM)
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));
	rowindices = new int[fm];
	colindices = new int[fn];
	for (int i = 0; i < fm; ++i)
		rowindices[i] = i;
	for (int i = 0; i < fn; ++i)
		colindices[i] = i;
	classname = "FCmatrixFull";
}


FCmatrixFull::FCmatrixFull(std::vector<Vec> v) : FCmatrix(v) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!(&v))
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));
	int fm = v.size();
	int fn = v[0].size();
	colindices = new int[fn];
	rowindices = new int[fm];
	for (int i = 0; i < fm; ++i)
		rowindices[i] = i;
	for (int i = 0; i < fn; ++i)
		colindices[i] = i;
	classname = "FCmatrixFull";
}

FCmatrixFull::FCmatrixFull(const FCmatrix& matrix){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this == &matrix)
	{
		std::cout << "It seems you know the wei..." << std::endl;
		return;
	}else if(!(&matrix))
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));


	int fm = matrix.GetRowN();
	int fn = matrix.GetColN();
	colindices = new int[fn];
	rowindices = new int[fm];
	for (int i = 0; i < fm; ++i)
		rowindices[i] = i;
	for (int i = 0; i < fn; ++i)
		colindices[i] = i;
	classname = "FCmatrixFull";
	M3.clear();
	M3.resize(fm, Vec());
	for(int i=0; i < fm; ++i)
 		M3[i] = matrix[i];
}

FCmatrixFull::FCmatrixFull(int fm, int fn, double a){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	classname = "FCmatrixFull";
	colindices = new int[fn];
	rowindices = new int[fm];
	for (int i = 0; i < fm; ++i)
		rowindices[i] = i;
	for (int i = 0; i < fn; ++i)
		colindices[i] = i;
	for (int i = 0; i < fm; i++)
		M3.push_back(Vec(fn, a));
}


/*
FCmatrixFull::~FCmatrixFull(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	delete rowindices;
	delete colindices;
}
*/

//------------------------------------------
//-----------------METHODS------------------
//------------------------------------------

int FCmatrixFull::GetRowN() const{
	return M3.size();
}

int FCmatrixFull::GetColN() const{
	return M3[0].size();
}

Vec FCmatrixFull::GetRow(int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int size = this->GetRowN();
	double x[size];
	for (int j = 0; j < size; ++j)
	{
		x[j] = M3[rowindices[i]][j];
	}
	return Vec(size, x);
}

Vec FCmatrixFull::GetCol(int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int size = M3.size();
	double nums[size];
	for (int j = 0; j < size; ++j)
		nums[j] = M3[rowindices[j]][i];
	return Vec(size, nums);
}

int FCmatrixFull::GetRowMax(int i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	double max = fabs(M3[i][0]);
	int i_max = 0;
	for (int j = 0; j < M3[i].size(); ++j)
	{
		if(fabs(M3[i][j]) > max)
		{
			max = fabs(M3[i][j]);
			i_max = j;
		}
	}
	return i_max;
}

int FCmatrixFull::GetColMax(int j){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	Vec v = this->GetCol(j);
	double max = fabs(v[0]);
	int i_max = 0;
	for (int i = 0; i < v.size(); ++i)
	{
		if(fabs(v[i]) > max)
		{
			max = fabs(v[i]);
			i_max = i;
		}
	}
	return i_max;
}

double FCmatrixFull::Determinant(){
	int n = this->GetColN();
	if(n != this->GetRowN()) return 0;

	if(n == 2) return ( (*this)[0][0] * (*this)[1][1] -
						(*this)[0][1] * (*this)[1][0] );
	
	if(n == 3) return ( (*this)[0][0] * (*this)[1][1] * (*this)[2][2] +
						(*this)[0][1] * (*this)[1][2] * (*this)[2][0] +
						(*this)[1][0] * (*this)[2][1] * (*this)[0][2] -
						(*this)[0][2] * (*this)[1][1] * (*this)[2][0] -
						(*this)[0][1] * (*this)[1][0] * (*this)[2][2] -
						(*this)[0][0] * (*this)[1][2] * (*this)[2][1] );

	FCmatrixFull Dummy(*this);

	double determinant = FCmatrixAlgorithms::GaussEliminationSimple(Dummy);
	for(int i=0; i<n ; i++)
		determinant *= Dummy[i][i];

	return determinant;
}

void FCmatrixFull::swapRows(int a, int b){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int dumb = rowindices[a];
	rowindices[a] = rowindices[b];
	rowindices[b] = dumb;
}

void FCmatrixFull::swapCols(int a, int b){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	  for (int i = 0; i < this->GetRowN(); ++i)
	  		M3[i].swap(a,b);
	  int dumb = colindices[a];
	  colindices[a] = colindices[b];
	  colindices[b] = dumb;
}

std::vector<int> FCmatrixFull::GetRowIndices() const{
	return std::vector<int> (rowindices, rowindices + this->GetColN());
}

std::vector<int> FCmatrixFull::GetColIndices() const{
	return std::vector<int> (colindices, colindices + this->GetRowN());
}

double FCmatrixFull::Norm1(){
	double sum = 0;
	double dum = 0;
	for (int i = 0; i < this->GetColN(); i++)
	{
		dum = 0;
		for (int j = 0; j < this->GetRowN(); j++)
		{
			dum += fabs(M3[j][i]);
		}
		if(dum > sum)
			sum = dum;
	}
	return sum;
}

double FCmatrixFull::NormInf(){
	double sum = 0;
	double dum = 0;
	for (int i = 0; i < this->GetRowN(); i++)
	{
		dum = 0;
		for (int j = 0; j < this->GetColN(); j++)
		{
			dum += fabs(M3[i][j]);
		}
		if(dum > sum)
			sum = dum;
	}
	return sum;
}

//right now it's only for square matrixes
FCmatrixFull FCmatrixFull::GetTranspose(){
    FCmatrixFull solution(*this);

    //Filling solution-matrix
    for(int i = 0; i < (*this).GetRowN(); i++) {
        for(int j = 0; j <(*this)[0].size(); j++) {
            solution[j][i] = (*this)[i][j];
        }
    }
    
    return solution;
}


FCmatrixFull FCmatrixFull::GetInverse(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	int n = this->GetColN();
	if(n != this->GetRowN()){
		throw std::invalid_argument(Form("[%s] Can't invert a non square matrix!\n",  __PRETTY_FUNCTION__));
	}

	std::vector <Vec> X;
	Vec x;

	for(int i=0; i<n; i++){
		Vec b(n);
		b[i] = 1;

		EqSolver Sol(*this,b);
		x = Sol.GaussEliminationSolver();
		X.push_back(x);
	}
	FCmatrixFull Inverse(X);

	return Inverse.GetTranspose();
}


//--------------------------------------------
//-----------------OPERATORS------------------
//--------------------------------------------

//Copy Assignment
void FCmatrixFull::operator=(const FCmatrixFull& matrix) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(this == &matrix) return;

	int m = matrix.GetRowN();
	int n = matrix.GetColN();

	if(this->GetRowN() != matrix.GetRowN()){
		if(rowindices)
			delete rowindices;

		rowindices = new int[m];
	}
	if(this->GetColN() != matrix.GetColN()){
		if(colindices)
			delete colindices;

		colindices = new int[n];
	}
	if(m == 0 || n == 0) return;

	std::vector<int> dumb = matrix.GetRowIndices();
	for (int i = 0; i < m; ++i) rowindices[i] = dumb[i];
	dumb.clear();
	dumb = matrix.GetColIndices();
	for (int i = 0; i < n; ++i) colindices[i] = dumb[i];
	if(this->GetRowN() != m) M3.resize(m, Vec(n));
	for (int i = 0; i < m; ++i)
		M3[i] = matrix.M3[i];
}

void FCmatrixFull::operator=(const FCmatrixSparse& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!(&mtx))
		return;
	if(rowindices) delete rowindices;
	if(colindices) delete colindices;


	M3.clear();
	M3.resize(0);
	int m = mtx.GetRowN();
	int n = mtx.GetColN();

	//setting up fresh rowindices and colindices
	rowindices = new int[m];
	colindices = new int[n];

	for (int i = 0; i < m; ++i)
		rowindices[i] = i;
	for (int i = 0; i < n; ++i)
		colindices[i] = i;

	for (int i = 0; i < m; ++i)
		M3.push_back(mtx[i]);

}

Vec& FCmatrixFull::operator[](int i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(i < 0 || i >= M3.size())
		throw std::invalid_argument(Form("[%s] That number (%d) cant be reached!\n",  __PRETTY_FUNCTION__, i));

	return M3[rowindices[i]];
}

Vec FCmatrixFull::operator[](int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(i < 0 || i >= M3.size())
		throw std::invalid_argument(Form("[%s] That number (%d) cant be reached!\n",  __PRETTY_FUNCTION__, i));

	return M3[rowindices[i]];
}

FCmatrixFull FCmatrixFull::operator+(const FCmatrixFull& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] The matrices are not of equal sizes (%ld x %d) vs. (%d x %d)!\n",  __PRETTY_FUNCTION__, this->M3.size(), this->M3[0].size(), mtx.GetSize(), mtx[0].size()));

	std::vector<Vec> v;
	for (int i = 0; i < this->GetColN(); ++i)
		v.push_back(this->GetRow(i) + mtx[i]);
	return FCmatrixFull(v);
}

FCmatrixFull FCmatrixFull::operator-(const FCmatrixFull& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] The matrices are not of equal sizes (%ld x %d) vs. (%d x %d)!\n",  __PRETTY_FUNCTION__, this->M3.size(), this->M3[0].size(), mtx.GetSize(), mtx[0].size()));

	std::vector<Vec> v;
	for (int i = 0; i < this->GetColN(); ++i)
		v.push_back(this->GetRow(i) - mtx[i]);
	return FCmatrixFull(v);
}

void FCmatrixFull::operator-= (const FCmatrixFull& mtx) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] The matrices are not of equal sizes (%ld x %d) vs. (%d x %d)!\n",  __PRETTY_FUNCTION__, this->M3.size(), this->M3[0].size(), mtx.GetSize(), mtx[0].size()));
	for (int i = 0; i < this->GetColN(); ++i)
		M3[i] -= mtx[i];
}

void FCmatrixFull::operator+= (const FCmatrixFull& mtx) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] The matrices are not of equal sizes (%ld x %d) vs. (%d x %d)!\n",  __PRETTY_FUNCTION__, this->M3.size(), this->M3[0].size(), mtx.GetSize(), mtx[0].size()));
	for (int i = 0; i < this->GetColN(); ++i)
		M3[i] += mtx[i];
}

FCmatrixFull FCmatrixFull::operator*(const FCmatrixFull& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] You cant multiply these two matrices, the number of columns of the first (%d) should equal the number of rows of the second (%d)!", __PRETTY_FUNCTION__, this->GetColN(), mtx.GetRowN()));
	
	std::vector<Vec> v;
	int m = mtx.GetRowN();
	int n = mtx.GetColN();
	std::vector<Vec> v_coiso;
	double x[n];

	//arranging the FCmatrix in vectors, that can be used for dot products
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			x[j] = mtx[j][i];
		v_coiso.push_back(Vec(m, x));
	}


	//dot products
	for (int i = 0; i < this->GetRowN(); ++i)
	{
		double x[n];
		for (int j = 0; j < n; ++j)
			x[j] = this->GetRow(i).dot(v_coiso[j]);
		v.push_back(Vec(n, x));
	}


	return FCmatrixFull(v);
}

FCmatrixFull FCmatrixFull::operator*(double lambda){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!this)
		throw std::invalid_argument(Form("[%s]: Nullptr!", __PRETTY_FUNCTION__));

	std::vector<Vec> v;
	for (int i = 0; i < this->GetRowN(); ++i)
		v.push_back(lambda*(this->GetRow(i)));
	return FCmatrixFull(v);
}

void FCmatrixFull::operator*=(double lambda){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	for (int i = 0; i < this->GetRowN(); ++i)
		M3[i] = M3[i]*lambda;
}

void FCmatrixFull::operator*=(const FCmatrixFull& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] You cant multiply these two matrices, the number of columns of the first (%d) should equal the number of rows of the second (%d)!", __PRETTY_FUNCTION__, this->GetColN(), mtx.GetRowN()));
	
	std::vector<Vec> v;
	int m = mtx.GetRowN();
	int n = mtx.GetColN();
	double * x = new double[m];
	std::vector<Vec> v_coiso;

	//arranging the FCmatrix in vectors, that can be used for dot products
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			x[j] = mtx[j][i];
		v_coiso.push_back(Vec(m, x));
	}

	//dot products
	for (int i = 0; i < this->GetRowN(); ++i)
	{
		x = new double[n];
		for (int j = 0; j < n; ++j)
			x[j] = this->GetRow(i).dot(v_coiso[j]);
		v.push_back(Vec(n, x));
	}

	delete x;
	*this = FCmatrixFull(v);
}

Vec FCmatrixFull::operator*(const Vec& vec){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this->GetColN() != vec.size())
		throw std::invalid_argument(Form("[%s] You cant multiply these two, the number of columns of the matrix (%d) should equal the number of rows of the vector (%d)!", __PRETTY_FUNCTION__, this->GetColN(), vec.size()));
	double x[this->GetRowN()];
	for (int i = 0; i < this->GetRowN(); ++i)
		x[i] = this->GetRow(i).dot(vec);
	return Vec(this->GetRowN(), x);
}


//--------------------------------------------
//----------------PRINT METHODS--------------_
//--------------------------------------------

void FCmatrixFull::PrintIndices() const{
	std::cout << "Indices are as follows: " << std::endl;
	for (int i = 0; i < this->GetRowN(); ++i)
		std::cout << rowindices[i] << ", " << std::endl;
}

void FCmatrixFull::Print() const{
	std::cout << "MatrixFull: " << std::endl;
	for (int i = 0; i < this->GetRowN(); ++i)
		std::cout << "              " << (*this)[i] << std::endl;
}

std::ostream& operator<<(std::ostream& s, const FCmatrixFull& mtx){
	s << "Matrix Full: " << std::endl;
  	for (int i=0; i<mtx.GetRowN(); ++i)
  		s << "              " << mtx[i] << std::endl;

 	 return s;
}