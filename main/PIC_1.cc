#include <time.h>
#include <math.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH1.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TView.h"
#include "TPad.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TImage.h"

#include "constants.h"


using std::vector;

int main(int argc, char *argv[]) {
	clock_t tStart = clock();

	double L = 1;
	double A = 2;
	double omega = 1e12;
	double temp;

	//funcao densidade de carga
	TF1* ro_Q = new TF1("ro_Q", [&](double*x, double *p){ //x[0]-> posicao ; p[0]-> tempo
		return A*(cos(M_PI*(x[0]-cos(omega*p[0]))/(2*L)))*(cos(M_PI*(x[0]-cos(omega*p[0]))/(2*L)))/(M_PI+sin(M_PI*(cos(omega*p[0])+L/2)/(2*L))-sin(M_PI*(cos(omega*p[0])-L/2)/(2*L))); }
		, -L/2, L/2, 1);

	//funcao que é integrada p = integral(d*q)
	TF1* ro_Q_d = new TF1("ro_Q_d", [&](double*x, double *p){
		ro_Q->SetParameter(0,p[0]);
		return x[0]*ro_Q->Eval(x[0]); }, -L/2, L/2, 1);

	//funcao de momento dipolar eletrico
	TF1* elec_dip_mom = new TF1("elec_dip_mom", [&](double*x, double *p){
		ro_Q_d->SetParameter(0,x[0]); //tempo passa a ser o x (já não é funçao da posicao)
		return ro_Q_d->Integral(-L/2,L/2); }, 0, 100, 0);

	//griffiths 11.58
	TF1* E_theta = new TF1("E_theta", [&](double*x, double *p){ // x[0]-> theta ; p[0]-> t ; p[1]-> r
		return mu_zero*elec_dip_mom->Derivative2(p[0]-p[1]/c_light)*(sin(x[0])/p[1])/(4*M_PI); }, 0, 4*M_PI, 2);


	TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
	auto gr = new TGraph();

	E_theta->SetParameter(1,1); //raio 1 metro para os graficos
	E_theta->SetParameter(0,0); // tempo a 0

	/*
	for (int i=0; i<20; i++){
        c1->cd();
        E_theta->SetParameter(0,i*1e-14);
        for (int j=0; j<100; j++) gr->SetPoint(j,2*M_PI*j/100, E_theta->Eval(2*M_PI*j/100));
        gr->Draw();
    	cout << i << endl;
        c1->Print("simulation.gif+5");
    }
	c1->Print("simulation.gif++");
	*/

	for(int i=0; i<200; i++){
		E_theta->SetParameter(0,i*1e-14); // tempo a 0
		gr->SetPoint(i,i*1e-14,E_theta->Eval(M_PI/2));
	}

	cout << tStart - clock() << endl;

	gr->Draw();
	c1->SaveAs("Files_Images_PIC/E_t.pdf");

	return 0;
}