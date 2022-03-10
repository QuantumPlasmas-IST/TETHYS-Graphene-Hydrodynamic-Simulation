


#include "includes/Fluid1DLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/DyakonovShurBoundaryLib.h"

#include "TethysBaseLib.h"
#include "includes/TethysMathLib.h"

#include "includes/Cell1DLib.h"
#include "includes/StateVecLib.h"




#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;


int main(int argc, char **argv){


StateVec U1(2,3);
StateVec U2(4,5);
StateVec U(U1);
U=U2/3;
cout << "\n\nTEsting operators"<<endl;
	cout << U <<endl;
cout <<"\n\n";

SetUpParameters parameters(2, 1, 0, 0, 0, 0,0, 1, 1);

Fluid1D teste(parameters);

teste.CflCondition();
cout <<"S\t"<<teste.GetVelSnd()<<endl;
cout <<"VF\t"<<teste.GetVelFer()<<endl;
cout <<"dt\t"<<teste.GetDt()<<endl;


//teste.SetSound([](float x) { return 1.0f+0.5f*tanh(8.0f*cos(2.0f*5.0f*MAT_PI*x)); });
teste.SetSound();


///////////////////
//testing algorithm
teste.CreateFluidFile();
teste.CreateHdf5File();

teste.SaveSound();



teste.InitialCondTest();
//teste.InitialCondRand();

	for (float h = 0.0f; h < 2.5f ; h+=teste.GetDt()) {
		teste.WriteFluidFile(h);
		if (GrapheneFluid1D::TimeStepCounter % 3 == 0) {
			teste.CopyFields();
			teste.SaveSnapShot();
		}
		GrapheneFluid1D::TimeStepCounter++;


		//teste.Richtmyer();
		teste.RungeKuttaTVD();
		//teste.McCormack();





BoundaryCondition::XPeriodic(teste);
//DyakonovShurBoundaryCondition::DyakonovShurBc(teste);
//DirichletBoundaryCondition::Density(teste,1.0f,1.0f);
//DirichletBoundaryCondition::VelocityX(teste,.1f,0.1f);


	}
	teste.WriteAttributes();
	teste.CloseHdf5File();


	int Nx=151;
	StateVec *Utest;
	Utest = new StateVec[Nx]();
	for (int i = 0; i < Nx; i++ ){
		Utest[i].n()=1.0;
		Utest[i].v()=(i>Nx/3 && i<2*Nx/3 ) ? 1.0f : 0.1f;
	}


	ofstream outputfile;
	outputfile.open ("TVDtest.dat");
	for (int i = 0; i < Nx; ++i) {
		CellHandler1D cell(i, nullptr, Utest);
		outputfile <<i <<"\t"<<  Utest[i]  <<"\t"<< cell.VanLeer(Utest,i) <<endl;
	}

	outputfile.close();

///////////////////



/*
	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	GrapheneFluid2D graph(parameters);
*/
	/*......CFL routine to determine dt...............................*/
/*
	graph.CflCondition();
	graph.SetDt(graph.GetDt()*0.25f);
	dt=graph.GetDt();
*/
	/*................................................................*/

	/*.........Fixed or variable vel_snd value........................*/

	//graph.SetSound();
	//graph.SetSimulationTime();
	//graph.SetTmax(5.0f);

	/*................................................................*/

	/*.........Output files and streams...............................*/
	//graph.CreateFluidFile();
	//graph.CreateHdf5File();
	 /*
	if(parameters.SaveMode){
		graph.SaveSound();
	}*/
	/*................................................................*/



	//GrapheneFluid2D::BannerDisplay();
	//graph.WelcomeScreen();

	/*...............Initialization...................................*/
	//graph.InitialCondWave();
	//graph.InitialCondRand();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
//	BoundaryCondition::SetSlope(0.0f);
//	BoundaryCondition::SetBottomEdge(graph);
//	BoundaryCondition::SetTopEdge(graph);
	/*................................................................*/
/*
	int Nx=150,Ny=150;
	float dx=1.0f/Nx;
	float dy=1.0f/Ny;
	float * source = new float[Nx*Ny]();
	float * laplacian = new float[Nx*Ny]();
	float * divergence = new float[Nx*Ny]();
	float * gradient_x = new float[Nx*Ny]();
	float * gradient_y = new float[Nx*Ny]();

	for (int i = 0; i < Nx ; ++i) {
		for (int j = 0; j < Ny ; ++j) {
			int k = i+j*Nx;
			float x=i*dx;
			float y=j*dx;
			//source[k] = pow((x - 0.5)*(y - 0.5),2.0);
			source[k] =sin(4*x/(3*y + 0.5));
		}
	}


	MathUtils::LaplacianField(source,laplacian,dx,Nx,Ny);
	MathUtils::GradientField(source,gradient_x,gradient_y,dx,dy,Nx,Ny);

	std::string previewfile1 = "2D_laplacian.dat" ;
	std::ofstream data_preview1;
	data_preview1.open (previewfile1);
	data_preview1 << scientific;

	std::string previewfile2 = "2D_gradientX.dat" ;
	std::ofstream data_preview2;
	data_preview2.open (previewfile2);
	data_preview2 << scientific;

	std::string previewfile3 = "2D_gradientY.dat" ;
	std::ofstream data_preview3;
	data_preview3.open (previewfile3);
	data_preview3 << scientific;


	for (int i = 0; i < Nx ; ++i) {
		for (int j = 0; j < Ny ; ++j) {
			int k = i+j*Nx;
			data_preview1 << laplacian[k] <<"\t";
			data_preview2 << gradient_x[k] <<"\t";
			data_preview3 << gradient_y[k] <<"\t";
		}
		data_preview1<<"\n";
		data_preview2<<"\n";
		data_preview3<<"\n";
	}
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Nx=150;

/*
	gsl_matrix * FDmatrix2;
	FDmatrix2 = gsl_matrix_calloc (Nx, Nx) ;
	float kappa=1/(dx*dx);
	float coeff;
	for (int i=0;i<Nx;i++) {
		for (int j = 0; j < Nx; j++) {
			if(j==i){
				gsl_matrix_set(FDmatrix2, i, j, -2*kappa);
			}else
			if(j==i+1){
				gsl_matrix_set(FDmatrix2, i, j, 1*kappa);
			}else
			if(j==i-1){
				gsl_matrix_set(FDmatrix2, i, j, 1*kappa);
			}
		}
	}
	//primeira linha 2	−5	4	−1
	gsl_matrix_set(FDmatrix2, 0, 0, 2*kappa);
	gsl_matrix_set(FDmatrix2, 0, 1, -5*kappa);
	gsl_matrix_set(FDmatrix2, 0, 2, 4*kappa);
	gsl_matrix_set(FDmatrix2, 0, 3, -1*kappa);
	//ultima linha
	gsl_matrix_set(FDmatrix2, Nx-1, Nx-1, -2*kappa);
	gsl_matrix_set(FDmatrix2, Nx-1, Nx-1-1, 5*kappa);
	gsl_matrix_set(FDmatrix2, Nx-1, Nx-1-2, -4*kappa);
	gsl_matrix_set(FDmatrix2, Nx-1, Nx-1-3, 1*kappa);
*/
/*
	gsl_matrix * FDmatrix3 ;
	FDmatrix3 = gsl_matrix_calloc (Nx, Nx) ;
	float kappa=1/(dx*dx*dx);
	for (int i=0;i<Nx;i++) {
		for (int j = 0; j < Nx; j++) {
			if(j==i){
				gsl_matrix_set(FDmatrix3, i, j, 0.0*kappa);
			}else
			if(j==i+1){
				gsl_matrix_set(FDmatrix3, i, j, -1*kappa);
			}else
			if(j==i-1){
				gsl_matrix_set(FDmatrix3, i, j, 1*kappa);
			}else
			if(j==i-2){
				gsl_matrix_set(FDmatrix3, i, j, -0.5*kappa);
			}else
			if(j==i+2){
				gsl_matrix_set(FDmatrix3, i, j, 0.5*kappa);
			}
		}
	}
//primeira linha 2	−5	4	−1
	gsl_matrix_set(FDmatrix3, 0, 0, -2.5*kappa);
	gsl_matrix_set(FDmatrix3, 0, 1, 9*kappa);
	gsl_matrix_set(FDmatrix3, 0, 2, -12*kappa);
	gsl_matrix_set(FDmatrix3, 0, 3, 7*kappa);
	gsl_matrix_set(FDmatrix3, 0, 4, -1.5*kappa);
//segunda linha
	gsl_matrix_set(FDmatrix3, 1, 0, 0*kappa);
	gsl_matrix_set(FDmatrix3, 1, 1, -2.5*kappa);
	gsl_matrix_set(FDmatrix3, 1, 2, 9*kappa);
	gsl_matrix_set(FDmatrix3, 1, 3, -12*kappa);
	gsl_matrix_set(FDmatrix3, 1, 4, 7*kappa);
	gsl_matrix_set(FDmatrix3, 1, 5, -1.5*kappa);
//ultima linha
	gsl_matrix_set(FDmatrix3, Nx-1, Nx-1, 2.5*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1, Nx-1-1, -9*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1, Nx-1-2, 12*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1, Nx-1-3, -7*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1, Nx-1-4, 1.5*kappa);
//penultima linha
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1, 0*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1-1, 2.5*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1-2, -9*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1-3, 12*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1-4, -7*kappa);
	gsl_matrix_set(FDmatrix3, Nx-1-1, Nx-1-5, 1.5*kappa);
*/

/*
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nx; ++j) {
			cout <<  gsl_matrix_get(FDmatrix2,i,j) <<"\t";
		}
		cout << "\n";
	}
*/
/*
	gsl_vector *sol = gsl_vector_alloc (Nx);
	gsl_vector *rhs = gsl_vector_alloc (Nx);

	for (int i=0;i<Nx;i++){
		gsl_vector_set(rhs,i,0.0f+0.05/cosh(10.0*(i*dx-.5))); //copiar o array de floats para o vectordouble
	}
	gsl_blas_dgemv(CblasNoTrans, 1, FDmatrix3, rhs, 0.0, sol);

	FILE *fp;
	fp = fopen("3rdDerivative.txt", "w");
	gsl_vector_fprintf (fp, sol, "%g");

	//gsl_permutation_free (p);
	gsl_vector_free (sol);
*/
	return 0;
}