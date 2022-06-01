#include "FCmatrixSparse.h"

//#define DEBUG

//--------------------------------------------
//-----------------CONSTRUCTORS---------------
//--------------------------------------------

FCmatrixSparse::FCmatrixSparse(){
#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif	
	classname = "FCmatrixSparse";
}


//Warning for someone reading the following code:
//We setted up these constructors so it will intake 
//a matrix as it is, with all the zeros and stuff
//because we didnt find it normal for someone to do the
//counting all by itself. Have a good trip!

FCmatrixSparse::FCmatrixSparse(double ** fM, int fm, int fn) : n(fn) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
    if(!fM)
		throw std::invalid_argument("null pointer!\n");
	
	std::cout << "[FCmatrixSparse] Warning: You cant change the values once they are here!" << std::endl;

	//now this is the real deal
	double x[fm*fn];
	double bs[fm];
	double cs[fm*fn];
	int    ctr  = 0;
	int    bs_ctr = 0;
	bool   bs_flag = true;
	int    zero_ctr = 0;
	for (int i = 0; i < fm; ++i)
	{
		bs_flag = true;
		zero_ctr = 0;
		for (int j = 0; j < fn; ++j)
		{
			if(fM[i][j] != 0) 
			{
				x[ctr] = fM[i][j];
				cs[ctr] = j;
				if(bs_flag == true){
					bs_flag = false;
					bs[bs_ctr] = ctr;
					bs_ctr++;
				}
				ctr++;
			}
			else{
				if(zero_ctr == n-1)
				{
					bs[bs_ctr] = -1;
					bs_ctr++;
				}
				zero_ctr++;
			}
		}
	}
 

	//so.... this part of the code is here because of the following
	//if a matrix has more than one row full of zeros, and we're gonna
	//force it, in the bs array, to have the same value as the index right after that one
	//then if we have more than one row full of zeros "in a row", than we must do
	//a while loop, just to check if everything is setted up nice and easy
	bool bs_cond = true;

	while(bs_cond){

		bs_cond = false;
		for (int i = 0; i < fm; ++i)
		{
			if(bs[i] == -1 && i+1 < fm){
				bs[i] = bs[i+1];
				bs_cond = true;
			}
			else if(bs[i] == -1){
				bs[i] = ctr;
				bs_cond = true;
			}
		}
	}

	M3.push_back(Vec(ctr, x));
	M3.push_back(Vec(bs_ctr, bs));
	M3.push_back(Vec(ctr, cs));
	classname = "FCmatrixSparse";
}

FCmatrixSparse::FCmatrixSparse(double * fM_i, int fm, int fn) : n(fn) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
    if(!fM_i)
		throw std::invalid_argument("null pointer!\n");
	
	if(fm < 1 || fn < 1)
		throw std::invalid_argument("Non positive number of rows or columns!");

	std::cout << "[FCmatrixSparse] Warning: You cant change the values once they are here!" << std::endl;

	double fM[fm][fn];
	int    plus = 0;

	//setting up a similar thing to the previous constructor
	//in order to make things easy :D we physicists do it the easy way
	for (int i = 0; i < fm; ++i)
	{
		for (int j = 0; j < fn; ++j)
		{
			fM[i][j] = fM_i[plus];
			plus++;
		}
	}

	//now this is the real deal
	double x[fm*fn];
	double bs[fm];
	double cs[fm*fn];
	int    ctr  = 0;
	int    bs_ctr = 0;
	bool   bs_flag = true;
	int    zero_ctr = 0;
	for (int i = 0; i < fm; ++i)
	{
		bs_flag = true;
		zero_ctr = 0;
		for (int j = 0; j < fn; ++j)
		{
			if(fM[i][j] != 0) 
			{
				x[ctr] = fM[i][j];
				cs[ctr] = j;
				if(bs_flag == true){
					bs_flag = false;
					bs[bs_ctr] = ctr;
					bs_ctr++;
				}
				ctr++;
			}
			else{
				if(zero_ctr == n-1)
				{
					bs[bs_ctr] = -1;
					bs_ctr++;
				}
				zero_ctr++;
			}
		}
	}
 

	//so.... this part of the code is here because of the following
	//if a matrix has more than one row full of zeros, and we're gonna
	//force it, in the bs array, to have the same value as the index right after that one
	//then if we have more than one row full of zeros "in a row", than we must do
	//a while loop, just to check if everything is setted up nice and easy
	bool bs_cond = true;

	while(bs_cond){

		bs_cond = false;
		for (int i = 0; i < fm; ++i)
		{
			if(bs[i] == -1 && i+1 < fm){
				bs[i] = bs[i+1];
				bs_cond = true;
			}
			else if(bs[i] == -1){
				bs[i] = ctr;
				bs_cond = true;
			}
		}
	}

	M3.push_back(Vec(ctr, x));
	M3.push_back(Vec(bs_ctr, bs));
	M3.push_back(Vec(ctr, cs));
	classname = "FCmatrixSparse";
}

FCmatrixSparse::FCmatrixSparse(std::vector<Vec> v){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!&v)
		throw std::invalid_argument("null pointer!\n");

	std::cout << "[FCmatrixSparse] Warning: You cant change the values once they are here!" << std::endl;
	
	int fm = v.size();
	int fn = v[0].size();

	//now this is the real deal
	double x[fm*fn];
	double bs[fm];
	double cs[fm*fn];
	int    ctr  = 0;
	int    bs_ctr = 0;
	bool   bs_flag = true;
	int    zero_ctr = 0;
	for (int i = 0; i < fm; ++i)
	{
		bs_flag = true;
		zero_ctr = 0;
		for (int j = 0; j < fn; ++j)
		{
			if(v[i][j] != 0) 
			{
				x[ctr] = v[i][j];
				cs[ctr] = j;
				if(bs_flag == true){
					bs_flag = false;
					bs[bs_ctr] = ctr;
					bs_ctr++;
				}
				ctr++;
			}
			else{
				if(zero_ctr == fn-1)
				{
					bs[bs_ctr] = -1;
					bs_ctr++;
				}
				zero_ctr++;
			}
		}
	}
 

	//so.... this part of the code is here because of the following
	//if a matrix has more than one row full of zeros, and we're gonna
	//force it, in the bs array, to have the same value as the index right after that one
	//then if we have more than one row full of zeros "in a row", than we must do
	//a while loop, just to check if everything is setted up nice and easy
	bool bs_cond = true;

	while(bs_cond){

		bs_cond = false;
		for (int i = 0; i < fm; ++i)
		{
			if(bs[i] == -1 && i+1 < fm){
				bs[i] = bs[i+1];
				bs_cond = true;
			}
			else if(bs[i] == -1){
				bs[i] = ctr;
				bs_cond = true;
			}
		}
	}

	n = fn;
	M3.push_back(Vec(ctr, x));
	M3.push_back(Vec(bs_ctr, bs));
	M3.push_back(Vec(ctr, cs));
	classname = "FCmatrixSparse";
}

FCmatrixSparse::FCmatrixSparse(const FCmatrixSparse& mtx){
	n = mtx.GetColN();
	M3.clear();
	M3.resize(mtx.M3.size(), Vec());
	M3 = mtx.M3;
}

//--------------------------------------------
//-----------------METHODS--------------------
//--------------------------------------------

Vec FCmatrixSparse::GetA() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	//returns the vector with all the values of the !=0 entries
	return M3[0];
}

Vec FCmatrixSparse::GetB() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	//returns the vector with the indices of the first values != 0 in each row, relatively to A
	return M3[1];
}

Vec FCmatrixSparse::GetC() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	//returns the vector with all the column indices of the != values
	return M3[2];
}


Vec FCmatrixSparse::GetRow(int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	return (*this)[i];
}

Vec FCmatrixSparse::GetCol(int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	double x[M3[1].size()];
	for (int j = 0; j < M3[1].size(); ++j)
		x[j] = (*this)[j][i];
	return Vec(M3[1].size(),x);
}

int FCmatrixSparse::GetRowMax(int i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(i < 0 || i >= n)
		throw std::invalid_argument("Row must be non negative or smaller than the size of the matrix");
	int index = M3[1][i];
	double max = M3[0][index];
	int cols_size = 0;

	if(i+1 < M3[1].size()) cols_size = M3[1][i+1] - M3[1][i];
	else cols_size = M3[0].size() - M3[1][i];

	for (int j = 0; j < cols_size; ++j)
		if(M3[0][index+j] > max)
			max = M3[0][index+j];

	for (int j = 0; j < n; ++j)
		if((*this)[i][j] == max)
			return j;

	std::cout << "Seems like this thing aint working good!" << std::endl;
	return 0;
}

int FCmatrixSparse::GetColMax(int j){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(j < 0 || j >= n)
		throw std::invalid_argument("Column must be non negative or smaller than the size of the matrix!\n");

	int index = 0;
	double max = 0;
	int ctr = 0;

	//getting the maximum value
	for (int i = 0; i < M3[2].size(); ++i)
		if(M3[2][i] == j)
			if(M3[0][i] > max)
					max = M3[0][i];

	for (int i = 0; i < M3[1].size(); ++i)
		if((*this)[i][j] == max)
			return i;

	std::cout << "Seems like this thing aint working good!" << std::endl;
	return 0;
}

int  FCmatrixSparse::GetRowN() const{return M3[1].size();}
int  FCmatrixSparse::GetColN() const{return n;}

//these methods are abandoned in this type of matrix: it would be to heavy computationally
//and for that we have MatrixFull, which does this much better
std::vector<int> FCmatrixSparse::GetRowIndices() const{
	std::vector<int> v;
	//this method is in development, for now it is abandoned
	return v;
}

std::vector<int> FCmatrixSparse::GetColIndices() const{
	std::vector<int> v;
	//this method is in development, for now it is abandoned
	return v;
}

double FCmatrixSparse::Determinant(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	//this method is in development, for now it is abandoned
	return 0;
}

void FCmatrixSparse::SetA(const Vec& a){
	if(M3.size() == 0)
		M3.push_back(a);
	else
		M3[0] = a;
}

void FCmatrixSparse::SetB(const Vec& b){
	if(M3.size() == 1)
		M3.push_back(b);
	else
		M3[1] = b;
}

void FCmatrixSparse::SetC(const Vec& c){
	if(M3.size() == 2)
		M3.push_back(c);
	else
		M3[2] = c;
}

//--------------------------------------------
//-----------------OPERATORS------------------
//--------------------------------------------

void FCmatrixSparse::operator=(const FCmatrixSparse& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	M3.clear();
	M3.resize(mtx.M3.size(), Vec());
	M3[0] = mtx.M3[0];
	M3[1] = mtx.M3[1];
	M3[2] = mtx.M3[2];
	n = mtx.n;
}

Vec& FCmatrixSparse::operator[](int i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int ctr = M3[1][i];
	//creating an all 0 Vec
	Vec row(n);

	//this work for all lines except the last one
	if(i+1 < M3[1].size()){
		for (int j = M3[1][i]; j < M3[1][i+1]; ++j)
		{
			//this gives us the amount of non 0 elements in the row
			row[M3[2][ctr]] = M3[0][ctr];
			//std::cout << "Aight: " << row[M3[2][ctr]] << std::endl;
			ctr++;
		}
	}
	
	//the last line, at last, ah ah! the irony!
	//its size will be... hum... well
	//it will start in B[i]. And it shall last for well, A.size()-B[i].
	else{
		for (int j = M3[1][i]; j < M3[0].size(); ++j)
		{
			row[M3[2][ctr]] = M3[0][ctr];
			//std::cout << "Not so good: " << row[M3[2][ctr]] << std::endl;
			ctr++;
		}
	}
	Vec * b = new Vec(row);
	return *b;
}

Vec FCmatrixSparse::operator[](int i) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int ctr = M3[1][i];
	Vec row(n);
	if(i+1 < M3[1].size()){
		for (int j = M3[1][i]; j < M3[1][i+1]; ++j)
		{
			row[M3[2][ctr]] = M3[0][ctr];
			ctr++;
		}
	}
	else{
		for (int j = M3[1][i]; j < M3[0].size(); ++j)
		{
			row[M3[2][ctr]] = M3[0][ctr];
			ctr++;
		}
	}
	return row;
}

FCmatrixSparse FCmatrixSparse::operator*(double lambda){
	FCmatrixSparse a;
	a.SetA(this->GetA()*lambda);
	a.SetB(this->GetB());
	a.SetC(this->GetC());
	a.n = n;
	return a;
}

FCmatrixSparse FCmatrixSparse::operator+(const FCmatrix& mtx){
	if(!&mtx)
		throw std::invalid_argument("Nullpointer!!\n");
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument("Different size matrices!!\n");

	std::vector<Vec> a;
	for (int i = 0; i < this->GetRowN(); ++i)
		a.push_back((*this)[i] + mtx[i]);
	return FCmatrixSparse(a);
}

FCmatrixSparse FCmatrixSparse::operator-(const FCmatrix& mtx){
	if(!&mtx)
		throw std::invalid_argument("Nullpointer!!\n");
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument("Different size matrices!!\n");

	std::vector<Vec> a;
	for (int i = 0; i < this->GetRowN(); ++i)
		a.push_back((*this)[i] - mtx[i]);
	return FCmatrixSparse(a);
}

void FCmatrixSparse::operator-=(const FCmatrix& mtx){
	if(!&mtx)
		throw std::invalid_argument("Nullpointer!!\n");
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument("Different size matrices!!\n");
	std::vector<Vec> a;
	for (int i = 0; i < this->GetRowN(); ++i)
		a.push_back((*this)[i] - mtx[i]);
	(*this) = FCmatrixSparse(a);
}

void FCmatrixSparse::operator+=(const FCmatrix& mtx){
	if(!&mtx)
		throw std::invalid_argument("Nullpointer!!\n");
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument("Different size matrices!!\n");

	std::vector<Vec> a;
	for (int i = 0; i < this->GetRowN(); ++i)
		a.push_back((*this)[i] + mtx[i]);
	(*this) = FCmatrixSparse(a);
}


//--------------------------------------------
//-----------------PRINTS---------------------
//--------------------------------------------


void FCmatrixSparse::Print() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	std::vector<Vec> dumb;
	int ctr = 0;

	for (int i = 0; i < M3[1].size(); ++i)
	{
		//creating an all 0 Vec
		Vec row(n);

		//this work for all lines except the last one
		if(i+1 < M3[1].size()){
			for (int j = M3[1][i]; j < M3[1][i+1]; ++j)
			{
				//this gives us the amount of non 0 elements in the row
				row[M3[2][ctr]] = M3[0][ctr];
				ctr++;
			}
			dumb.push_back(row);
		}
		//the last line, at last, ah ah! the irony!
		//its size will be... hum... well
		//it will start in B[i]. And it shall last for well, A.size()-B[i].
		else{
			for (int j = M3[1][i]; j < M3[0].size(); ++j)
			{
				row[M3[2][ctr]] = M3[0][ctr];
				ctr++;
			}
			dumb.push_back(row);
		}
	}

	std::cout << "MatrixSparse:" << std::endl;
	for (int i = 0; i < dumb.size(); ++i)
		std::cout << "              " << dumb[i] << std::endl;
}

std::ostream& operator<<(std::ostream& s, const FCmatrixSparse& mtx){
	s << "MatrixSparse:" << std::endl;
	for (int i = 0; i < mtx.GetRowN(); ++i)
		s << "              " << mtx[i] << std::endl;
	return s;
}