#include "Vec.h"

#include "TROOT.h"

//#define DEBUG

//--------------------------------------------
//-----------------CONSTRUCTORS---------------
//--------------------------------------------

Vec::Vec(int i, double x) :  N(i) {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if (N<0) throw std::invalid_argument(Form("[%s] received negative number of elements...!\n", __PRETTY_FUNCTION__));
  entries = new double[N];
  indices = new int[N];
  for(int j=0; j<N; j++)
  	indices[j] = j;

  std::fill_n(entries, N, x);
}

Vec::Vec(int i, const double* x) : Vec(i, 0.) { //c++11 on...
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if (x)
    std::copy(x, x+i, entries);
  else  
    throw std::invalid_argument(Form("[%s] null pointer to array...!\n", __PRETTY_FUNCTION__));
}

Vec::Vec(int i, const int* x) : Vec(i, 0.) { //c++11 on...
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if (x)
    std::copy(x, x+i, entries);
  else  
    throw std::invalid_argument(Form("[%s] null pointer to array...!\n", __PRETTY_FUNCTION__));
}

Vec::Vec(const Vec& v) : Vec(v.N, v.entries) { //c++11 on...
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
}

Vec::Vec(std::vector<double> b): N(b.size()) {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  entries = new double[N];
  indices = new int[N];
  for(int j=0; j<N; j++)
  {
    indices[j] = j;
    entries[j] = b[j];
  }
}


//--------------------------------------------
//-----------------DESTRUCTOR-----------------
//--------------------------------------------

Vec::~Vec() {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif 
  delete entries;
  delete indices;
}


//--------------------------------------------
//-----------------OPERATORS------------------
//--------------------------------------------

//Read a member of Vec
double Vec::operator[](int i) const {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if (i>=N) 
    throw std::invalid_argument(Form("[%s] index out of bounds...(i=%d N=%d)!\n", __PRETTY_FUNCTION__, i, N));  
  return entries[indices[i]];
}

//Change a member of Vec
double& Vec::operator[](int i) {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if (i>=N) 
    throw std::invalid_argument(Form("[%s] index out of bounds...(i=%d N=%d)!\n", __PRETTY_FUNCTION__, i, N));  
  
  return entries[indices[i]];
}

//Copy assignment
void Vec::operator=(const Vec& v){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
	if (this != &v){
    if (v.N != N){
      N = v.N;
      delete entries;
      delete indices;
      entries = new double[N];
      indices = new int[N];
    }
    std::copy(v.entries, v.entries+N, entries);   
    for(int i=0; i<N; i++) 
      indices[i] = v.indices[i];
	}
}

//Adding 2 Vec's
Vec Vec::operator+(const Vec& other){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (other.size()!=N) 
    throw std::invalid_argument(Form("[%s] Vector size [%d] different from input vector size[%d]\n", __PRETTY_FUNCTION__, N, other.size()));  

  double* resultado = new double[N];
  for(int i=0; i<N; i++)
  	resultado[i] = entries[indices[i]] + other[i];

  Vec Res(N, resultado);
  delete resultado;

  return Res;
}

//Adding 2 Vec's, using +=
void Vec::operator+= (const Vec& v) { 
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (v.N != N)
    throw std::invalid_argument(Form("[%s] objects with different size...(N=%d v.N=%d)!\n", __PRETTY_FUNCTION__, N, v.N));        

  for (int i=0; i<N; ++i)
    entries[indices[i]] += v[i];
}

//Subtracting two Vec's
Vec Vec::operator-(const Vec& other){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (other.size()!=N) 
    throw std::invalid_argument(Form("[%s] Vector size [%d] different from input vector size[%d]\n", __PRETTY_FUNCTION__, N, other.size()));  

  double* resultado = new double[N];
  for(int i=0; i<N; i++)
  	resultado[i] = entries[indices[i]] - other[i];

  Vec Res(N, resultado);
  delete resultado;

  return Res;
}

Vec Vec::operator-(){
  double* resultado = new double[N];
  for(int i=0; i<N; i++)
    resultado[i] = -entries[indices[i]];

  Vec Res(N, resultado);
  delete resultado;

  return Res;
}

//Subtracting two Vec's, using -=
void Vec::operator-= (const Vec& v) { 
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (v.N != N) 
    throw std::invalid_argument(Form("[%s] objects with different size...(N=%d v.N=%d)!\n", __PRETTY_FUNCTION__, N, v.N));        
  
  for (int i=0; i<N; ++i) 
    entries[indices[i]] -= v[v.indices[i]];
}

//Multiplying Vec with scalar
Vec Vec::operator*(const double x) const {
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (abs(x-1.)<1E-9)
    return *this;

  double a[N];
  for (int i=0; i<N; ++i) 
    a[i] = entries[indices[i]] * x;

  return Vec(N, a);
}

//Multiplying each entry of a Vec by the corresponding entry in the other Vec
Vec Vec::operator*(const Vec& other) const{
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (other.size()!=N) 
    throw std::invalid_argument(Form("[%s] Vector size [%d] different from input vector size[%d]\n", __PRETTY_FUNCTION__, N, other.size()));  
  
  double* resultado = new double[N];
  for(int i=0; i<N; i++) 
    resultado[i] = entries[indices[i]] * other[i];

  Vec Res(N, resultado);
  delete resultado;

  return Res;
}

//Mult. but using *=
void Vec::operator*= (const Vec& v) { 
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (v.N != N)
    throw std::invalid_argument(Form("[%s] objects with different size...(N=%d v.N=%d)!\n", __PRETTY_FUNCTION__, N, v.N));
  
  for (int i=0; i<N; ++i) 
    entries[indices[i]] *= v[i];
}

//Returns the Normalized Vec
Vec Vec::operator!(){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  double resultado[N];
  double norma = this->mod();
  for(int i=0; i<N; i++) 
    resultado[i] = entries[indices[i]]/norma;

  return Vec(N, resultado);
}


//--------------------------------------------
//-----------------FRIEND METHODS-------------
//--------------------------------------------

//cout << Vec << endl;
std::ostream& operator<<(std::ostream& s, const Vec& v) {
  s << "[";
  for (int i=0; i<v.N; ++i) {
    s << std::fixed << std::setprecision(6) << v[i];
    if (i<v.N-1) s << ", ";
  }
  s << "]";

  return s;
}


Vec operator*(double x, const Vec& v){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (abs(x-1.)<1E-9)
    return Vec(v.size(), 0.);

  double a[v.size()];
  for (int i=0; i<v.size(); ++i) 
    a[i] = v[v.indices[i]] * x;

  return Vec(v.size(), a);
}


//------------------------------------------
//-----------------METHODS------------------
//------------------------------------------

//Dot product between 2 Vec's
double Vec::dot(const Vec& v){
	if(v.N != N)
		throw std::invalid_argument(Form("[%s] objects with different size (%d , %d)", __PRETTY_FUNCTION__, N, v.N));
	double dumb = 0;
  for (int i = 0; i < N; i++)
    dumb += (*this)[i]*v[i];
  
  return dumb;	
}

//External product between 2 Vec's
Vec Vec::ex(const Vec& v){
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif

  	//3D external product
	if(N==3 && v.N==3){
		double pe[3]={entries[indices[1]]*v[2] - entries[indices[2]]*v[1] ,  entries[indices[2]]*v[0] - entries[indices[0]]*v[2], entries[indices[0]]*v[1] - entries[indices[1]]*v[0]};
		return Vec(3, pe);
	}
  //7D external product
	else if(N==7 && v.N==7){
		double pe[7]= {
			entries[indices[1]]*v[3] - entries[indices[3]]*v[1] + entries[indices[2]]*v[6] - entries[indices[6]]*v[2] + entries[indices[4]]*v[5] - entries[indices[5]]*v[4],
			entries[indices[2]]*v[4] - entries[indices[4]]*v[2] + entries[indices[3]]*v[0] - entries[indices[0]]*v[3] + entries[indices[5]]*v[6] - entries[indices[6]]*v[5],
			entries[indices[3]]*v[5] - entries[indices[5]]*v[3] + entries[indices[4]]*v[1] - entries[indices[1]]*v[4] + entries[indices[6]]*v[0] - entries[indices[0]]*v[6],
			entries[indices[4]]*v[6] - entries[indices[6]]*v[4] + entries[indices[5]]*v[2] - entries[indices[2]]*v[5] + entries[indices[0]]*v[1] - entries[indices[1]]*v[0],
			entries[indices[5]]*v[0] - entries[indices[0]]*v[5] + entries[indices[6]]*v[3] - entries[indices[3]]*v[6] + entries[indices[1]]*v[2] - entries[indices[2]]*v[1],
			entries[indices[6]]*v[1] - entries[indices[1]]*v[6] + entries[indices[0]]*v[4] - entries[indices[4]]*v[0] + entries[indices[2]]*v[3] - entries[indices[3]]*v[2],
			entries[indices[0]]*v[2] - entries[indices[2]]*v[0] + entries[indices[1]]*v[5] - entries[indices[5]]*v[1] + entries[indices[3]]*v[4] - entries[indices[4]]*v[3],};
		return Vec(7, pe);
	}
	else
		throw std::invalid_argument(Form("[%s] x product not defined for dimensions diffent than 3 and 7. Your input was (%d,%d)", __PRETTY_FUNCTION__, N, v.N));
}

//swap two numbers (Actually, it only swaps the indices, but shhh)
void Vec::swap(int i, int j){
	if(std::max(i,j) >= N)
		throw std::invalid_argument(Form("[%s]indices out of range (%d , %d)", __PRETTY_FUNCTION__, N, std::max(i,j)));
	
  if (i!=j){
    int a = indices[j];
    indices[j] = indices[i];
    indices[i] = a;
  }
}

double Vec::sumAbs(const Vec& v){
	if(v.N != N)
		throw std::invalid_argument(Form("[%s] objects with different size (%d , %d)", __PRETTY_FUNCTION__, N, v.N));
	
  return std::accumulate(entries, entries+N, 0, 
		[](double accum, double x){return accum+fabs(x);});	
}

//Size of Vec
int Vec::size() const{ return this->N;}

//Returns the sqrt of the sum of squares of each member in Vec
double Vec::mod() const{
	double norma=0;
	for(int i=0; i<N; i++) 
    norma+=entries[indices[i]]*entries[indices[i]];

	return sqrt(norma);
}

double Vec::AbsMax() const{
  double dumb = fabs((*this)[0]);
  for (int i = 0; i < N; i++)
    if(fabs((*this)[i]) > dumb)
      dumb = fabs((*this)[i]);
  return dumb;
}

double* Vec::data(){
  return entries;
}

//Prints Vec
void Vec::Print() {
	std::cout << *this << std::endl;
}

//Setting Values in a Vec, after it has been declared
void Vec::SetEntries (int n, double* x){	
#ifdef DEBUG
  printf("[%s]\n", __PRETTY_FUNCTION__ );
#endif
  if (N==0){
  	delete indices;
  	indices = new int[n];
  	for(int i=0; i<n; i++) indices[i] = i;
  }
	if(n!= N){
		delete[] entries;
		entries = new double[n];
		N = n;
	}

	for (int i = 0; i < n; ++i) 
    (*this)[i] = x[i];
}

int* Vec::GetIndices(){
  return indices;
}
//--------------------------------------------
//-----------------FOREIGNERS-----------------
//--------------------------------------------

//Swapping 2 Vec's
void swap(Vec& v1, Vec& v2){
	if(v1.size() != v2.size())
		throw std::invalid_argument(Form("[%s]Vectors must be of same size(%d , %d)", __PRETTY_FUNCTION__, v1.size(), v2.size()));
	
  for(int i=0; i<v1.size(); i++) 
    std::swap(v1[i], v2[i]);
}