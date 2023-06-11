/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/



#include "includes/InitialConditionLib.h"
#include "includes/SetUpParametersLib.h"
#include "includes/DirichletBoundaryLib.h"
#include "includes/DyakonovShurBoundaryLib.h"
#include "includes/GrapheneFluid2DLib.h"
#include "TethysBaseLib.h"

#ifndef MAT_PI
#	define MAT_PI 3.14159265358979323846
#endif


using namespace std;

float ff_top(float x){
	return 0.125*x;
}

float ff_bottom(float x){
	return 0.125*x;
}

int main(int argc, char **argv){

	SetUpParameters parameters(argc, argv);
	parameters.DefineGeometry();

	float t=0.0;
	float dt;		// time step

	int NX = 401;
    int NY = 201;
	Geometry Geom(NX,NY);
	Geom.fronteira.D.set_Domain(ff_top,ff_bottom);
    Geom.fronteira.set_Edge();
	
	Geom.dominio.dom = Geom.fronteira.D.dom;
/*	cout << "evaluating" << endl;
	for(int i = 0; i < Geom.fronteira.D.dom.size(); i++){
		if(Geom.fronteira.edg[i] == true){
			if(Geom.fronteira.edg[i-1] == false && Geom.fronteira.edg[i+1] == false && Geom.fronteira.edg[i-NX] == false && Geom.fronteira.edg[i+NX] == false ){
				cout << "wait wait wait" << endl;
			}
			//cout << "edgint[" << i << "] = " << Geom.fronteira.edgint[i] << "      "; 
		}
	}*/

//	cout << "edgint size  " << Geom.fronteira.edgint.size() << endl;

/*	for(int i = 0; i < Geom.fronteira.edg.size(); i++){
		cout << "edg[" << i << "] = " << Geom.fronteira.edg[i] << "      "; 
	}
*/
	GrapheneFluid2D graph(parameters);

	/*......CFL routine to determine dt...............................*/
	graph.CflCondition();
	//graph.SetDt(graph.GetDt()*0.25f);
	dt=graph.GetDt();
	/*................................................................*/
	
	/*.........Fixed or variable vel_snd value........................*/
	float sound=graph.GetVelSnd();
	float ly=graph.GetLengthY();
	//graph.GetLengthX();

	//std::function<float(float,float)> variationS = [=](float x,float y){ return sound+5.f* tanh(6.0f*cos(2.0f*MAT_PI*2.0f*x)); };
	//std::function<float(float,float)> variationS = [=](float x,float y){ return sound*(1+0.3f*x* abs(y-0.5*ly)); };
	//graph.SetSound(variationS);
	graph.SetSound();
	//graph.SetSimulationTime();
	//graph.SetTmax(3.f);

	/*................................................................*/

	/*.........Output files and streams...............................*/
	graph.CreateFluidFile();
	graph.CreateHdf5File();
	if(parameters.SaveMode){
		graph.SaveSound();
	}
	/*................................................................*/



	GrapheneFluid2D::BannerDisplay();
	graph.WelcomeScreen();

	/*...............Initialization...................................*/
//	graph.InitialCondGeneral([](float x,float y) { return 1.0f+0.1f/cosh(10.0f*sqrt((x-.5f)*(x-.5f)+(y-.5f)*(y-.5f))); },[](float x,float y) { return 0.5f; },[](float x,float y) { return 0.0f; });
//	graph.InitialCondGeneral([](float x,float y) { return 1.0f+0.1f/cosh(10.0f*(x-.5f)); },[](float x,float y) { return 0.5f/cosh(10.0f*(x-.5f)); },[](float x,float y) { return 0.0f; });
//	graph.InitialCondGeneral([](float x,float y) { return 0.8; },[](float x,float y) { return 0.5f; },[](float x,float y) { return 0.0f; });

	//graph.InitialCondRand();

//InitialCondition::Rand(graph);
//InitialCondition::InitialCondUniform(graph);
InitialCondition::InitialCondGeneral(graph, [](float x,float y) {return 1.0f;},[](float x,float y) { return 0.005f; },[](float x,float y) { return 0.0f; });

/*for(int k = 0; k < 201*401; k++){
	cout << graph.Umain[k].px() << endl;
}*/
//	graph.InitialCondTest();
	/*................................................................*/

	/*................Setting.the.lateral.boundaries..................*/
	//BoundaryCondition::SetSlope(0.0f);
	//BoundaryCondition::SetBottomEdge(graph);
	//BoundaryCondition::SetTopEdge(graph);
	/*................................................................*/
	
//	DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
	DirichletBoundaryCondition::XFreeLeft(graph);
	DirichletBoundaryCondition::XFreeRight(graph, &Geom);
//	DirichletBoundaryCondition::YClosedFreeSlip(graph);
	DirichletBoundaryCondition::YClosedNoSlipG(graph, &Geom);
//	BoundaryCondition::XPeriodic(graph);
//	BoundaryCondition::YPeriodic(graph);

	cout << "\033[1;7;5;33m Program Running \033[0m"<<endl;
	while (t <= graph.GetTmax()){
		int percentage=100*GrapheneFluid2D::TimeStepCounter/(graph.GetTmax()/dt);
		cout << percentage<<"%\033[?25l" << endl; //prints the percentage of simulation completed
//		cout << " t = " << t << endl;
		t += dt;
		GrapheneFluid2D::TimeStepCounter++;
		graph.Richtmyer(&Geom);
/*		if(t >= 0 && t <= 0.000003){
			bool debug_flag = 0;
			float fails[201*401];
			cout << "debugging" << endl;
			for(int k = 0; k < 201*401; k++){
				//cout << graph.Umain[k].px() << endl;
				if(!isfinite(graph.Umain[k].px()) || !isfinite(graph.Umain[k].n())){
					cout << "oops - Umain[" << k << "] = " << graph.Umain[k].n() << "      ";
					cout << k << "       ";
					debug_flag = 1;
				}
			}
			cout << endl;
			if(debug_flag == 1){
				cout << "oops"<< endl;
			}else{
				cout << "seems okay" << endl;
			}
		}*/
//		graph.Richtmyer();
		/*+++++++++++++++++++++++++++++++++++++*
		 * Change the boundary conditions here *
		 *+++++++++++++++++++++++++++++++++++++*/
//		DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
		DirichletBoundaryCondition::XFreeLeft(graph);
		DirichletBoundaryCondition::XFreeRight(graph, &Geom);
		DirichletBoundaryCondition::YClosedNoSlipG(graph,&Geom);
		//DirichletBoundaryCondition::YClosedFreeSlip(graph);

//		BoundaryCondition::XPeriodic(graph);
//		BoundaryCondition::YPeriodic(graph);


	//	if(graph.GetThermDiff()!=0.0){
	//		DirichletBoundaryCondition::Temperature(graph,0.22f, 0.22f, 0.22f, 0.22f);  // 300K corresponds to 0.22*Fermi temperature
	//	}
		if(graph.GetKinVis()!=0.0f ) {
			graph.ParabolicOperatorWeightedExplicit19('V',&Geom);
			//*+++++++++++++++++++++++++++++++++++++*
			// * Change the boundary conditions here *
			// *+++++++++++++++++++++++++++++++++++++
			//DyakonovShurBoundaryCondition::DyakonovShurBc(graph);
			DirichletBoundaryCondition::XFreeLeft(graph);
			DirichletBoundaryCondition::XFreeRight(graph, &Geom);
			DirichletBoundaryCondition::YClosedNoSlipG(graph, &Geom);
			//DirichletBoundaryCondition::YClosedFreeSlip(graph);
		}

/*		if(t >= 0.000185714){
			cout << "values debugging" << endl;
			for (int k = 0; k<NX*NY-NX-NY; k++){
				//cout << graph.Umain[k] << " ";
				if(!isfinite(graph.Umain[k].px())){
					for(int j = 0; j < 3; j++){
						cout << "Umain[" << k - NX - 1 + j << "].px = " << graph.Umain[k/NX - NX - 1 + j].px() << endl;
						cout << "Umain[" << k - 1 + j<< "].px = " << graph.Umain[k/NX - 1 + j].px() << endl;
						cout << "Umain[" << k + NX - 1 + j<< "].px = " << graph.Umain[k/NX + NX - 1 + j].px() << endl;  	
					}
				}
			}
		}*/

		//Record full hdf5 data
		if (parameters.SaveMode  && graph.Snapshot()) {
			graph.SaveSnapShot();
		}

		if( !( GrapheneFluid2D::TimeStepCounter % 2) ){
			graph.WriteFluidFile(t);
		}
		cout <<"\033[1G\033[2K"; //clears percentage of completion
	}
	//Record atributes on hdf5 file
	if(parameters.SaveMode) {
		graph.WriteAttributes();
	}
	graph.CloseHdf5File();


	cout << "\033[1A\033[2K\033[1;32mDONE!\033[0m\n";
	cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}