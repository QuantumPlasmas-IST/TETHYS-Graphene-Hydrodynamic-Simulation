#ifndef LaxWen_h
#define LaxWen_h

#include<iostream>

using namespace std;

/**********************************************************
*  This Class does not take into account E conservation   *
*                                                         *
*  3 conservative variables: density (den);               *
*                            momentumX (px= den*vX)       *
*                            momentumY (py= den*vY)       *
*                                                         *
**********************************************************/

/*******************************************************************************************
*        		In function SetBV(unsigned short, float*) : 							   *				
*																						   *
* int value specifies along which direction and for which variable's BV is being specified *
*																						   *
*		First number - boundary - 0->x0 | 1 -> x1 | 2 -> y0 | 3 -> y1 					   *
*		Second number - which variable - 0->den | 1 -> px | 2 -> py 					   *
*																						   *
*******************************************************************************************/

/******* Boundary Conditions *********
*  Indexes: bv[i] and bv_mat[i]      *
* 									 *
*  i = 0 -> x=0, den 				 *
*  i = 1 -> x=0, px					 *
*  ...								 *
*  i = 5 -> x=Nx*dx, py 		 	 *
*  i = 6 -> y=0, den 				 *
*  ... 								 *
*  i=11 -> y=Ny*dy, py				 *
* 									 *
* 									 *
*  bv_mat[i] size : Ny (i<6)		 *
*					Nx (i>5)		 *
* 									 *
* 									 *
*  bv[i] = false -> Dirichelet       *
* 		   true -> free end          *	
* 									 *
**************************************/

/********* BC in corners of grid - Hierarchy ************
*														*
*	Neumann+Neumann -> Adjacent corner 					*
*	Neumann+Dirichelet -> Dirichelet 					*
*	Dirichelet(0)+Dirichelet -> Dirichelet(0)			*
*	Dirichelet+Dirichelet -> Average					*
*														*
********************************************************/

/********+********* FLux Functions *******************
*													 * 
*	var_in and var_out are arrays of 3 floats		 *
*													 * 
*	var_in : variables at [i][j]					 *
*	var_out : flux at [i][j]						 *
*													 * 
*****************************************************/

/********** MacCormack ***********
*						         *
* in void MacCor(..., bool)      *
*						         *
* bool = true -> fix mcx & mcy   *
* 		 false -> rotation       *
*						         *
*********************************/

/********** var_cur ************
*							   *
*   var_cur[i][j][k]           *
* 							   *
*   i -> x cell (0 - Nx-1)	   *
*	j -> y cell	(0 - Ny-1)	   *
*	k -> variable (0,1,2)	   *
* 							   *
*******************************/

class LaxWen2D{
	public:
		LaxWen2D(); // Default constructor (1x1 grid, 100 steps, dt=.001, Dirichelet Boundary Conditions everywhere)
		LaxWen2D(float, int, float, int, float, bool* = NULL);
		~LaxWen2D(); // Destructor 
		void SetGrid(float, float, float, float);
		void SetTimeStep(float);
		void SetBV(bool*);
		void SetBV(unsigned short, float* var_bv); 
		void Rich2(float*** var_cur, void FluxX(float* var_in, float* var_out), void FluxY(float* var_in, float* var_out)); // Richtmeyer 2-step
		void SetMCbool(bool, bool);
		void MacCor(float*** var_cur, void FluxX(float* var_in, float* var_out), void FluxY(float* var_in, float *var_out), bool); // MacCormack Predictor-Corrector	
		void DumpVariables();
		void DumpBV(unsigned short);								 
	private:
		float dx, dy; // x and y step
		int Nx, Ny; // L = Nx*dx   W = Ny*dy
		float dt; // time step
		bool* bv; // Boundary conditions - false -> Dirichelet   |   true -> Von Neumann 0
		float** bv_mat;
		bool mcx, mcy; // false - predictor with i, i-1 | true - predictor step with i+1, i
};

#endif