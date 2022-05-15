#include "FCmatrix.h"

//#define DEBUG

FCmatrix::FCmatrix() {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
}

FCmatrix::FCmatrix(double** a, int m_s, int n_s){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!a)
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));

	if(m_s < 1 || n_s < 1)
		throw std::invalid_argument(Form("[%s] Invalid number of columns or rows!\n",  __PRETTY_FUNCTION__));
	//setting up M3, vector<Vec>
	for(int i=0; i<m_s; ++i)
		M3.emplace_back(n_s, a[i]);
}

FCmatrix::FCmatrix(double* a, int m_s, int n_s){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!a)
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));
	if(m_s < 1 || n_s < 1)
		throw std::invalid_argument(Form("[%s] Invalid number of columns or rows!\n",  __PRETTY_FUNCTION__));
	for(int i=0; i<m_s; ++i)
		M3.emplace_back(n_s, &a[i*n_s]);
}

FCmatrix::FCmatrix(const std::vector<Vec>& v){
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(!(&v))
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));

	for(int i=0; i< v.size(); ++i)
		M3.emplace_back(v[i]);
}

FCmatrix::FCmatrix(const FCmatrix& matrix) {
	#ifdef DEBUG
	  printf("[%s]\n", __PRETTY_FUNCTION__ );
	#endif
	if(this == &matrix)
	{
		std::cout << "It seems you know the wei..." << std::endl;
		return;
	}else if(!(&matrix))
		throw std::invalid_argument(Form("[%s] null pointer!\n",  __PRETTY_FUNCTION__));

	for(int i=0; i < matrix.GetRowN(); ++i)
	 	M3.emplace_back(matrix.M3[i]);
}

/////////////// methods

int FCmatrix::GetSize() const{
	return M3.size();
}
std::string FCmatrix::GetClassName() const{
	return classname;
}

void FCmatrix::Print() const{
	for (int i = 0; i < M3.size(); ++i)
	{
		for (int j = 0; j < M3[0].size(); ++j)
		{
			std::cout << M3[i][j] << " " << std::endl;
		}
		std::cout << "\n" << std::endl;
	}
}

/*
Vec& FCmatrix::operator[](int i){
	if(i < 0 || i >= M3.size())
		throw std::invalid_argument(Form("[%s] That number (%d) cant be reached!\n",  __PRETTY_FUNCTION__, i));
	return M3[i];
}

Vec FCmatrix::operator[](int i) const{
	if(i < 0 || i >= M3.size())
		throw std::invalid_argument(Form("[%s] That number (%d) cant be reached!\n",  __PRETTY_FUNCTION__, i));
	return M3[i];
} */
