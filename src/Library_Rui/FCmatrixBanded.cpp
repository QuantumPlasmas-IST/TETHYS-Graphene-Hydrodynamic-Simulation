#include "FCmatrixBanded.h"

//#define DEBUG


//-----------------------------------------------
//-----------------CONSTRUCTORS------------------
//-----------------------------------------------

FCmatrixBanded::FCmatrixBanded(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	classname = "FCmatrixBanded";
}

FCmatrixBanded::FCmatrixBanded (double ** fBands, int flower, int fupper, int fnDiag) : N(fnDiag), upperD(fupper), lowerD(flower) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	//we assume that fBands is an array with the bands of the matrix only
	classname = "FCmatrixBanded";
	if(N < 1)
		throw std::invalid_argument(Form("[%s] Amount of numbers in the main diagonal must be non negative or number of diagonals must be greater than 0 (%d)!!\n",  __PRETTY_FUNCTION__,fnDiag));

	int index = 0;
	for (int i = 0; i <= upperD; ++i)
	{
		M3.push_back(Vec(N-upperD+i, fBands[index]));
		index++;
	}
	for (int i = 0; i < lowerD; ++i)
	{
		M3.push_back(Vec(N-1-i, fBands[index]));
		index++;
	}
}

FCmatrixBanded::FCmatrixBanded(double * fBands, int flower, int fupper, int fnDiag) : N(fnDiag), upperD(fupper), lowerD(flower) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(fnDiag < 1 || flower < 0 || fupper < 0)
		throw std::invalid_argument(Form("[%s] Amount of numbers in the main diagonal must be non negative or number of diagonals must be greater than 0 (%d)!!\n",  __PRETTY_FUNCTION__,fnDiag));
   
    classname = "FCmatrixBanded";
   	int length = 0;
   	int index = 0;
   	//filling up the upper diagonals
   	for (int i = 0; i <= upperD; ++i)
   	{
   		length = N-upperD+i;
   		double dumb[length];
   		for (int k = 0; k < length; ++k) dumb[k] = fBands[index+k];
   		index += length;
   		M3.push_back(Vec(length, dumb));
   	}

   	for (int i = 0; i < lowerD; ++i)
   	{
   		length = N - 1 - i;
   		double dumb[length];
   		for (int k = 0; k < length; ++k) dumb[k] = fBands[index+k];
   		index += length;
   		M3.push_back(Vec(length, dumb));
   	}
}

FCmatrixBanded::FCmatrixBanded(std::vector<Vec>& bands){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	//checking the main diagonal
    int index = 0;
    for(int i= 0; i < bands.size(); ++i)
    	if(bands[index].size() < bands[i].size())
    		index = i;

    classname = "FCmatrixBanded";
    M3 = bands;
    N = bands[index].size();
    upperD = index;
    lowerD = bands.size() - upperD - 1;
}

FCmatrixBanded::FCmatrixBanded(const FCmatrixBanded& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(&mtx == this)
		return;
	lowerD = mtx.GetLowerD();
	upperD = mtx.GetUpperD();
	N = mtx.GetColN();
	M3.clear();
	M3.resize(mtx.M3.size(), Vec());
	M3 = mtx.M3;
}

FCmatrixBanded::FCmatrixBanded(const FCmatrix& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(&mtx == this)
		return;
	if(mtx.GetColN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] The matrix you want to change to matrix is not square!!\n",  __PRETTY_FUNCTION__));
	
	//in development
	upperD = 0;
	lowerD = 0;
	N = this->GetRowN();
	for (int i = 1; i < this->GetRowN(); i++)
	{
		if((*this)[0][i] != 0)
			upperD++;

		if((*this)[i][0] != 0)
			lowerD++;		
	}
	

}

//-----------------------------------------------
//-----------------METHODS-----------------------
//-----------------------------------------------

double FCmatrixBanded::Determinant(){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	return 0;
}

int FCmatrixBanded::GetRowN() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif	
	return N;
}

int FCmatrixBanded::GetColN() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	return N;
}

int FCmatrixBanded::GetRowMax(int i){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	Vec dum = this->GetRow(i);
	int max_index = 0;
	for (int j = 0; j < dum.size(); ++j)
		if(fabs(dum[j]) > fabs(dum[max_index]))
			max_index = j;
	return max_index;
}

int FCmatrixBanded::GetColMax(int j){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	int max_index = 0;
	double max = 0;
	for (int i = 0; i < N; ++i)
		if(fabs(this->GetRow(i)[j]) > max)
		{
			max = fabs(this->GetRow(i)[j]);
			max_index = i;
		}

	return max_index;
}

int FCmatrixBanded::GetUpperD()const{
	return upperD;
}

int FCmatrixBanded::GetLowerD()const{
	return lowerD;
}

Vec FCmatrixBanded::GetCol(int j) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	double x[N];
	for (int i = 0; i < N; ++i)
		x[i] = this->GetRow(i)[j];
	Vec * a = new Vec(N, x);
	return *a;
}

Vec FCmatrixBanded::GetRow(int j) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(j < 0 || j >= N)
  		throw std::invalid_argument(Form("[%s] Number of diagonals must be greater than 0 (%d)!!\n",  __PRETTY_FUNCTION__,N));
  
	double row[N];

	//filling the upper part of the matrix
	for (int i = 0; i < N - j; ++i)
	{	
		if(i > upperD) row[j+i] = 0;
		else if(upperD-i >= 0 && j < M3[upperD-i].size()) row[j+i] = M3[upperD-i][j];
		else row[j+i] = 0;
	}

	//filling the lower part of the matrix
	int ctr = 0;
	for (int i = j; i > 0; --i)
	{
		if(i > lowerD) row[ctr] = 0;
		else row[ctr] = M3[upperD+i][ctr];
		ctr++;
	}
	return Vec(N, row);
}

Vec FCmatrixBanded::GetBand(int j) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	return M3[j];
}

std::vector<int> FCmatrixBanded::GetRowIndices() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	std::vector<int> a;
	return a;
}

std::vector<int> FCmatrixBanded::GetColIndices() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	std::vector<int> a;
	return a;
}

int FCmatrixBanded::GetMainDiagonalIndex() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
 	int index = 0;
    for(int i= 0; i < this->M3.size(); ++i)
    	if(this->M3[index].size() < this->M3[i].size())
    		index = i;
    return index;
}

//-----------------------------------------------
//-----------------OPERATORS---------------------
//-----------------------------------------------

Vec FCmatrixBanded::operator[](int j) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	return this->GetRow(j);
}

Vec& FCmatrixBanded::operator[](int j){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(j < 0 || j >= N)
  		throw std::invalid_argument(Form("[%s] Number of diagonals must be greater than 0 (%d)!!\n",  __PRETTY_FUNCTION__,N));
  
	double row[N];
	//filling the upper part of the matrix
	for (int i = 0; i < N - j; ++i)
	{	
		if(i > upperD) row[j+i] = 0;
		else row[j+i] = M3[upperD-i][j];
	}

	//filling the lower part of the matrix
	int ctr = 0;
	for (int i = j; i > 0; --i)
	{
		if(i > lowerD) row[ctr] = 0;
		else row[ctr] = M3[upperD+i][ctr];
		ctr++;
	}

	Vec * b = new Vec(N,row);
	return *b;
}

void FCmatrixBanded::operator=(const FCmatrixBanded& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(&mtx == this)
		return;
	this->M3.clear();
	this->M3.resize(mtx.M3.size(), Vec());
	this->M3 = mtx.M3;
	N = mtx.GetRowN();
	upperD = mtx.GetUpperD();
	lowerD = mtx.GetLowerD();
}

FCmatrixBanded FCmatrixBanded::operator+(const FCmatrixBanded& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] Matrices with different size!!\n",  __PRETTY_FUNCTION__));
	
	int this_ctr = this->GetMainDiagonalIndex();
	int mtx_ctr = mtx.GetMainDiagonalIndex();
	std::vector<Vec> a;

	int min_upperd_index = (mtx.GetUpperD() >= this->GetUpperD()) ? this->GetUpperD() : mtx.GetUpperD();
	int min_lowerd_index = (mtx.GetLowerD() >= this->GetLowerD()) ? this->GetLowerD() : mtx.GetLowerD();

	//in case they have a different number of upper diagonals
	if(mtx.GetUpperD() > min_upperd_index)
		for (int i = 0; i < mtx.GetUpperD() - min_upperd_index; ++i)
			a.push_back(mtx.GetBand(i));
	else if(this->GetUpperD() > min_upperd_index)
		for (int i = 0; i < this->GetUpperD() - min_upperd_index; ++i)
			a.push_back(this->GetBand(i));

	//adding the above the main diagonal bands
	for (int i = 0; i <= min_upperd_index; ++i)
		a.push_back(mtx.GetBand(mtx_ctr+i-min_upperd_index) + (*this).GetBand(this_ctr+i-min_upperd_index));

	for (int i = 1; i <= min_lowerd_index; ++i)
		a.push_back(mtx.GetBand(mtx_ctr+i) + (*this).GetBand(this_ctr+i));

	if(mtx.GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= mtx.GetLowerD() - min_lowerd_index; ++i)
			a.push_back(mtx.GetBand(mtx_ctr+i));
	else if(this->GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= this->GetLowerD() - min_lowerd_index; ++i)
			a.push_back(this->GetBand(this_ctr+i));

	return FCmatrixBanded(a);
}

FCmatrixBanded FCmatrixBanded::operator-(const FCmatrixBanded& mtx){
	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] Matrices with different size!!\n",  __PRETTY_FUNCTION__));
	
	int this_ctr = this->GetMainDiagonalIndex();
	int mtx_ctr = mtx.GetMainDiagonalIndex();
	std::vector<Vec> a;

	int min_upperd_index = (mtx.GetUpperD() >= this->GetUpperD()) ? this->GetUpperD() : mtx.GetUpperD();
	int min_lowerd_index = (mtx.GetLowerD() >= this->GetLowerD()) ? this->GetLowerD() : mtx.GetLowerD();

	//in case they have a different number of upper diagonals
	if(mtx.GetUpperD() > min_upperd_index)
		for (int i = 0; i < mtx.GetUpperD() - min_upperd_index; ++i)
			a.push_back(-mtx.GetBand(i));
	else if(this->GetUpperD() > min_upperd_index)
		for (int i = 0; i < this->GetUpperD() - min_upperd_index; ++i)
			a.push_back(this->GetBand(i));

	//adding the above the main diagonal bands
	for (int i = 0; i <= min_upperd_index; ++i)
		a.push_back(-mtx.GetBand(mtx_ctr+i-min_upperd_index) + (*this).GetBand(this_ctr+i-min_upperd_index));

	for (int i = 1; i <= min_lowerd_index; ++i)
		a.push_back(-mtx.GetBand(mtx_ctr+i) + (*this).GetBand(this_ctr+i));


	if(mtx.GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= mtx.GetLowerD() - min_lowerd_index; ++i)
			a.push_back(-mtx.GetBand(mtx_ctr+i));
	else if(this->GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= this->GetLowerD() - min_lowerd_index; ++i)
			a.push_back(this->GetBand(this_ctr+i));

	return FCmatrixBanded(a);
}

void FCmatrixBanded::operator+=(const FCmatrixBanded& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] Matrices with different size!!\n",  __PRETTY_FUNCTION__));
	
	int this_ctr = this->GetMainDiagonalIndex();
	int mtx_ctr = mtx.GetMainDiagonalIndex();
	std::vector<Vec> a;

	int min_upperd_index = (mtx.GetUpperD() >= this->GetUpperD()) ? this->GetUpperD() : mtx.GetUpperD();
	int min_lowerd_index = (mtx.GetLowerD() >= this->GetLowerD()) ? this->GetLowerD() : mtx.GetLowerD();

	//in case they have a different number of upper diagonals
	if(mtx.GetUpperD() > min_upperd_index)
		for (int i = 0; i < mtx.GetUpperD() - min_upperd_index; ++i)
			a.push_back(mtx.GetBand(i));
	else if(this->GetUpperD() > min_upperd_index)
		for (int i = 0; i < this->GetUpperD() - min_upperd_index; ++i)
			a.push_back(this->GetBand(i));

	//adding the above the main diagonal bands
	for (int i = 0; i <= min_upperd_index; ++i)
		a.push_back(mtx.GetBand(mtx_ctr+i-min_upperd_index) + (*this).GetBand(this_ctr+i-min_upperd_index));

	for (int i = 1; i <= min_lowerd_index; ++i)
		a.push_back(mtx.GetBand(mtx_ctr+i) + (*this).GetBand(this_ctr+i));

	if(mtx.GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= mtx.GetLowerD() - min_lowerd_index; ++i)
			a.push_back(mtx.GetBand(mtx_ctr+i));
	else if(this->GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= this->GetLowerD() - min_lowerd_index; ++i)
			a.push_back(this->GetBand(this_ctr+i));

	(*this).M3.clear();
	(*this).M3.resize(a.size(), Vec());
	(*this).M3 = a;
}

void FCmatrixBanded::operator-=(const FCmatrixBanded& mtx){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif

	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));
	if(this->GetColN() != mtx.GetColN() || this->GetRowN() != mtx.GetRowN())
		throw std::invalid_argument(Form("[%s] Matrices with different size!!\n",  __PRETTY_FUNCTION__));
	
	int this_ctr = this->GetMainDiagonalIndex();
	int mtx_ctr = mtx.GetMainDiagonalIndex();
	std::vector<Vec> a;

	int min_upperd_index = (mtx.GetUpperD() >= this->GetUpperD()) ? this->GetUpperD() : mtx.GetUpperD();
	int min_lowerd_index = (mtx.GetLowerD() >= this->GetLowerD()) ? this->GetLowerD() : mtx.GetLowerD();

	//in case they have a different number of upper diagonals
	if(mtx.GetUpperD() > min_upperd_index)
		for (int i = 0; i < mtx.GetUpperD() - min_upperd_index; ++i)
			a.push_back(-mtx.GetBand(i));
	else if(this->GetUpperD() > min_upperd_index)
		for (int i = 0; i < this->GetUpperD() - min_upperd_index; ++i)
			a.push_back(this->GetBand(i));

	//adding the above the main diagonal bands
	for (int i = 0; i <= min_upperd_index; ++i)
		a.push_back(mtx.GetBand(mtx_ctr+i-min_upperd_index) + (*this).GetBand(this_ctr+i-min_upperd_index));

	for (int i = 1; i <= min_lowerd_index; ++i)
		a.push_back(-mtx.GetBand(mtx_ctr+i) + (*this).GetBand(this_ctr+i));

	if(mtx.GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= mtx.GetLowerD() - min_lowerd_index; ++i)
			a.push_back(-mtx.GetBand(mtx_ctr+i));
	else if(this->GetLowerD() > min_lowerd_index)
		for (int i = 1; i <= this->GetLowerD() - min_lowerd_index; ++i)
			a.push_back(this->GetBand(this_ctr+i));

	(*this).M3.clear();
	(*this).M3.resize(a.size(), Vec());
	(*this).M3 = a;
}

FCmatrixBanded FCmatrixBanded::operator*(const FCmatrix& mtx){
	if(!&mtx)
		throw std::invalid_argument(Form("[%s] Nullpointer!!\n",  __PRETTY_FUNCTION__));

	FCmatrixBanded a;
	return a;
}

FCmatrixBanded FCmatrixBanded::operator*(double lambda){
	std::vector<Vec> a;

	//copying the bands
	for (int i = 0; i < this->M3.size(); ++i)
		a.push_back(M3[i]*lambda);
	return FCmatrixBanded(a);
}

//-----------------------------------------------
//-----------------PRINTS------------------------
//-----------------------------------------------

void FCmatrixBanded::Print() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	std::cout << "MatrixBanded:" << std::endl;
	for (int i = 0; i < N; ++i)
		std::cout << "              " << this->GetRow(i) << std::endl;	
}

void FCmatrixBanded::PrintBands() const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	std::cout << "MatrixBanded->Bands:" << std::endl;
	for (int i = 0; i < M3.size(); ++i)
		std::cout << "              " << M3[i] << std::endl;
}

std::ostream& operator<<(std::ostream& s, const FCmatrixBanded& mtx){
	s << "MatrixBanded:" << std::endl;
	for (int i = 0; i < mtx.GetRowN(); ++i)
		s << "              " << mtx[i] << std::endl;
	return s;
}

//---------------------------------------------
//-----------------Solver----------------------
//---------------------------------------------

Vec FCmatrixBanded::Thomas(Vec d) const{
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	if(this->GetUpperD() != 1 || this->GetLowerD() != 1)
		throw std::invalid_argument(Form("[%s] Matrix is not Tridiagonal [#Above D=%d][#Bellow D=%d]\n", __PRETTY_FUNCTION__, this->GetUpperD(), this->GetLowerD()));

	if(d.size() != N)
		throw std::invalid_argument(Form("[%s] b vector doesn't have the same size as M [N=%d & size.b=%d]\n", __PRETTY_FUNCTION__, N, d.size()));


	int n = N-1;
	
	//creating Vec for bands vec
	Vec band_c(this->GetBand(0));
	Vec band_b(this->GetBand(1));
	Vec band_a(this->GetBand(2));
	Vec c_star(n);
	Vec d_star(d);

	if(fabs(band_b[0])==0)
		throw std::invalid_argument(Form("[%s] Cannot use this method because of division by 0 in the first element of diagonal\n", __PRETTY_FUNCTION__));

	c_star[0] = band_c[0] / band_b[0];
  	d_star[0] = d[0] / band_b[0];

  	double m;
  	for (int i=1; i<n; i++) {
	    m = band_b[i] - band_a[i-1] * c_star[i-1];

	    if(fabs(m)==0)
			throw std::invalid_argument(Form("[%s] Cannot use this method because of division by 0 in m[%f]\n", __PRETTY_FUNCTION__,m));

	    c_star[i] = band_c[i] / m;
	    d_star[i] = (d[i] - band_a[i-1] * d_star[i-1]) / m;
  	}

  	d_star[n] = (d_star[n] - band_a[n-1]*d_star[n-1]) /
  				 (band_b[n] - band_a[n-1]*c_star[n-1]);

    for (int i=n-1; i>=0; i--)
    	d_star[i] -= c_star[i] * d_star[i+1];

    return d_star;
}