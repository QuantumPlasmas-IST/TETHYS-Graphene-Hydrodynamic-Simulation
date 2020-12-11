/* Calling DGELS using row-major order */

#include <stdio.h>
#include <openblas/lapacke.h>
#include <openblas/cblas.h>
int main (int argc, const char * argv[])
{
	double a[5*3] = {1,2,1,2,3,4,3,5,2,4,2,5,5,4,3};
	double b[5*2] = {-10,-3,12,14,14,12,16,16,18,16};
	lapack_int info,m,n,lda,ldb,nrhs;
	int i,j;

	m = 5;
	n = 3;
	nrhs = 2;
	lda = 3;
	ldb = 2;

	info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,a,lda,b,ldb);

	for(i=0;i<n;i++)
	{
		for(j=0;j<nrhs;j++)
		{
			printf("%lf ",b[i+j*n]);
		}
		printf("\n");
	}
	return(info);
}