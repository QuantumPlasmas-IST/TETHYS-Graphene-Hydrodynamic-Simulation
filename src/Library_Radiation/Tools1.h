#ifndef TOOLS1_H
#define TOOLS1_H

#include <string>
#include <cstdio>
#include <math.h>

bool check_int (std::string in, double min = 0, double max = 0);
bool check_nat (std::string in, double min = 0, double max = 0);
bool check_double(std::string in, double min = 0, double max = 0);

long int factorial(long int n);

int nCr(int n, int r);

void del(int** ptr, int n);
void del(int*** ptr, int n, int m);

void del(float** ptr, int n);
void del(float*** ptr, int n, int m);

void del(double** ptr, int n);
void del(double*** ptr, int n, int m);

void sort2(int, double*, double*);

double GetMax(int, double*);
double GetMin(int, double*);

//sorting

void swap(int* a, int* b);
void swap(double* a, double* b);

int partition(int x[], int low, int high);
int partition(double x[], int low, int high);

void quicksort(int x[], int low, int high);
void quicksort(double x[], int low, int high);

int binarySearch(int x[], int low, int high, int element);
int binarySearch(double x[], int low, int high, double element);

#endif