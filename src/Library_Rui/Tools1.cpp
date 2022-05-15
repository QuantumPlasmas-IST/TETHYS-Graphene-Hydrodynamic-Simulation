#include "Tools1.h"
#include <iostream>

//checking if int
bool check_int (std::string in, double min, double max){
	int nspace=0;
	for (int i = 0; i < in.length(); i++){
   		if (isdigit(in[i]) == false){
   			if(in[i]==' ') nspace++;
	   		if(i!=nspace  && in[i]==45 || in[i]!=' ' && in[i]!=45){
				std::cout << "[TOOLS1]: This is not an int!" << std::endl;
				return false;
			}
      	}
    }
    //test range
    //option for no range check
    if(min == 0 && max == 0) return true;
    //check range
    int douplelganger = atoi(in.c_str());
    if( douplelganger < min || douplelganger > max){
    	std::cout << "[TOOLS1]: The int is not in the specified range!" << std::endl;
    	return false;
    }
	return true;
}

//checking if natural number
bool check_nat (std::string in, double min, double max){
	for (int i = 0; i < in.length(); i++){
   		if (isdigit(in[i]) == false){
   			//if(in[i]==' ') nspace++;
   			if(in[i]!=' '){
				std::cout << "[TOOLS1]: This is not a natural number!" << std::endl;
				return false;
			}
      	}
    }

    //test range
    //option for no range check
    if(min == 0 && max == 0) return true;
    //check range
    unsigned int douplelganger = atoi(in.c_str());
    if(douplelganger < min || douplelganger > max){
    	std::cout << "[TOOLS1]: The number is not in the specified range!" << std::endl;
    	return douplelganger;
    }

	return true;
}

//checking if double
bool check_double (std::string in, double min, double max){
	int nspace=0;
	int pv=0;
	for (int i = 0; i < in.length(); i++){
   		if (isdigit(in[i]) == false){
   			if(in[i]==' ') nspace++;
   			if(in[i] == 44 || in[i] == 46) pv++; //comma counter
   			if(in[i] !=13 && in[i] !=' '  && in[i] != 32 && in[i] != 44 && in[i] != 45 && in[i] != 46 || pv==2 || i!=nspace && in[i]==45){
				//std::cout << "what?" << (int) in[i] << " " << i << std::endl;
				std::cout << "[TOOLS1]: This is not a float/double!" << std::endl;
				return false;
      		}
      	}
    }

    //test range
    //option for no range check
    if(min == 0 && max == 0) return true;
    //check range
    double douplelganger = atof(in.c_str());
    if(douplelganger < min || douplelganger > max){
    	std::cout << "[TOOLS1]: The number is not in the specified range!" << std::endl;
    	return false;
    }

	return true;
}

//deleting pointers functions
void del(int** ptr, int n){
	if (n>0){
		for(int i=0; i<n ; i++){
			delete[] ptr[i];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: int** successfully deleted!" << std::endl;
	}
}

void del(int*** ptr, int n, int m){
	if (n>0 && m>0){
		for(int j=0; j<n; j++){
			for(int i=0; i<m ; i++){
				delete[] ptr[j][i];
			}
			delete[] ptr[j];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: int*** successfully deleted!" << std::endl;
	}
}

void del(float** ptr, int n){
	if (n>0){
		for(int i=0; i<n ; i++){
			delete[] ptr[i];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: float** successfully deleted!" << std::endl;
	}
}

void del(float*** ptr, int n, int m){
	if (n>0 && m>0){
		for(int j=0; j<n; j++){
			for(int i=0; i<m ; i++){
				delete[] ptr[j][i];
			}
			delete[] ptr[j];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: float*** successfully deleted!" << std::endl;
	}
}

void del(double** ptr, int n){
	if (n>0){
		for(int i=0; i<n ; i++){
			delete[] ptr[i];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: double** successfully deleted!" << std::endl;
	}
}

void del(double*** ptr, int n, int m){
	if (n>0 && m>0){
		for(int i=0; i<n; i++){
			for(int j=0; j<m ; j++){
				delete[] ptr[i][j];
			}
			delete[] ptr[i];
		}
		delete[] ptr;
		std::cout << "[TOOLS1]: double*** successfully deleted!" << std::endl;
	}
}


double GetMax(int n, double* a){

	double dumb = a[0];
	for (int i = 0; i < n; i++)
		if(a[i] > dumb)
			dumb = a[i];
	
	return dumb;
}


double GetMin(int n, double* a){

	double dumb = a[0];
	for (int i = 0; i < n; i++)
		if(a[i] < dumb)
			dumb = a[i];
	
	return dumb;
}

long int factorial(long int n){
	if(n < 1){
		//std::cout << "[TOOLS2] No can do" << std::endl;
		return 1;
	}

	if(n<30){
	long int dumb = 1;
	for (int i = 0; i < n; i++)
		dumb *= (n-i);

	return dumb;
	}else
		return sqrt(2*M_PI*n)*pow((double) n, (double) n)*(1+1/(12*n)+1/(288*n*n));
}

int nCr(int n, int r){
    return factorial(n) / (factorial(r) * factorial(n - r)); 
}

////////////  quicksort  ////////////
void swap(int* a, int* b){ 
    int t = *a; 
    *a = *b; 
    *b = t; 
}

void swap(double* a, double* b){ 
    double t = *a; 
    *a = *b; 
    *b = t; 
} 

int partition(int x[], int low, int high){
	int pivot = x[high];
	int i = low - 1;

	for(int j=low; j<=(high-1); j++){
		if(x[j] <= pivot){
			i++;
			swap(&x[i],&x[j]);
		}
	}
	swap(&x[i+1],&x[high]);

	return (i+1);
}

int partition(double x[], int low, int high){
	double pivot = x[high];
	int i = low - 1;

	for(int j=low; j<=(high-1); j++){
		if(x[j] <= pivot){
			i++;
			swap(&x[i],&x[j]);
		}
	}
	swap(&x[i+1],&x[high]);

	return (i+1);
}

void quicksort(int x[], int low, int high){
	if(low<high){
		int pi = partition(x,low,high);

		quicksort(x,low,pi-1);
		quicksort(x,pi+1,high);
	}
}

void quicksort(double x[], int low, int high){
	if(low<high){
		int pi = partition(x,low,high);

		quicksort(x,low,pi-1);
		quicksort(x,pi+1,high);
	}
}

//////////  binarySearch  //////////
int binarySearch(int x[], int low, int high, int element){
	while(low <= high){
		int m = low + (high - low)/2;

		if(x[m] == element) return m;
		if(x[m] < element) low = m + 1;
		else high = m - 1;
	}
	std::cout << "The element " << element << " could not be found" << std::endl;
	return -1;
}

int binarySearch(double x[], int low, int high, double element){
	while(low <= high){
		int m = low + (high - low)/2;

		//condition
		if(fabs(x[m] - element) < 0.5) return m;
		if(x[m] < element) low = m + 1;
		else high = m - 1;
	}
	std::cout << "The element " << element << " could not be found" << std::endl;
	return -1;
}
//////////////////////////////////