


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

/*
StateVec U1(2,3);
StateVec U2(4,5);
StateVec U(U1);
U=U2/3;
cout << "\n\nTEsting operators"<<endl;
	cout << U <<endl;
cout <<"\n\n";
*/




SetUpParameters parameters(13, 10, 0, 0.0, 0, 0,0, 1, 1);

GrapheneFluid1D teste(parameters);

teste.CflCondition();
teste.SetDt(teste.GetDt()*0.05);
cout <<"S\t"<<teste.GetVelSnd()<<endl;
cout <<"VF\t"<<teste.GetVelFer()<<endl;
cout <<"dt\t"<<teste.GetDt()<<endl;

teste.SetSound();


///////////////////
//testing algorithm
teste.CreateFluidFile();
teste.CreateHdf5File();

teste.SaveSound();



teste.InitialCondTest();
//teste.InitialCondRand();

	for (float h = 0.0f; h < 0.5f ; h+=teste.GetDt()) {
	//for (float h = 0.0f; h < 5*teste.GetDt() ; h+=teste.GetDt()) {
		teste.WriteFluidFile(h);
		if (GrapheneFluid1D::TimeStepCounter % 25 == 0) {
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
	outputfile.open ("TVD_limiter_test.dat");
	for (int i = 0; i < Nx; ++i) {
		CellHandler1D cell(i, nullptr, Utest);
		outputfile <<i <<"\t"<<  Utest[i]  <<"\t"<< cell.VanLeer(Utest,i) <<endl;
	}

	outputfile.close();


	return 0;
}