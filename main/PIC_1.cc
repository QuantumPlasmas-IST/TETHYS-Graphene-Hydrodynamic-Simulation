#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

//Não se chama TETHYS_2D_RadiationPattern.cpp ainda porque ainda não vai buscar os ficheiros input automaticamente (e porque é mais rapido compilar com um nome pequeno)
//Le de um ficheiro fixo, terá de um dia ler de um ficheiro criado numa sessão
//Usa ROOT apenas para visualização gráfico, mais tarde terá tudo a ser guardado em ficheiros e uma versão em que não precisa de root para correr
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphPolar.h>
#include <TGraphPolargram.h>
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

double c_light = 2.99792458e8;
double mu_zero = 4*M_PI*1e-7;
double epsilon_zero = 8.62607015e-34;

using namespace std;

int main(int argc, char *argv[]) {
	clock_t tStart = clock();

	string namefile = "Files_Images_PIC/electro_2D_S=21.00vF=10.50vis=0.010odd=0.000l=0.001wc=0.00therm=0.00.dat";

	ifstream infile1(namefile, ios::in);

    string line;
    double temp_d;
    int Nl = 112;
    int Nc = 15;

    vector<double> time;
    vector<double> edpX;
    vector<double> edpX_dd;
    vector<double> edpY;
    vector<double> edpY_dd;

    for(int i=0; i<Nl; i++){
        getline(infile1, line);
        stringstream ss;
        ss << line;
        string temp_s="";

        for(int j=0; j<Nc; j++){
            ss >> temp_s;
            if(j==0){
            	if (stringstream(temp_s) >> temp_d);
            	time.push_back(temp_d);
            }
            if(j==9){
            	if (stringstream(temp_s) >> temp_d);
            	edpX.push_back(temp_d);
            }
            if(j==11){
            	if (stringstream(temp_s) >> temp_d);
            	edpX_dd.push_back(temp_d);
            }
            if(j==12){
            	if (stringstream(temp_s) >> temp_d);
            	edpY.push_back(temp_d);
            }
            if(j==14){
            	if (stringstream(temp_s) >> temp_d);
            	edpY_dd.push_back(temp_d);
            }
        }
    }

    infile1.close();
    cout << "[FIM DA LEITURA DO FICHEIRO " << namefile << " ]" << endl;


    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
	auto gr_edpX = new TGraph();
	auto gr_edpY = new TGraph();

    for(int i=0; i<Nl; i++){
		gr_edpX->SetPoint(i,time[i],edpX[i]);
		gr_edpY->SetPoint(i,time[i],edpY[i]);
	}

	c1->cd();
	gr_edpX->Draw();
	c1->SaveAs("Files_Images_PIC/edpX_t.pdf");
	c1->Clear();
	gr_edpY->Draw();
	c1->SaveAs("Files_Images_PIC/edpY_t.pdf");


	auto Poyting = [&](int i, double x, double y, double z){ //modulo do vetor Poyting // tempo está no i
		double r2 = x*x+y*y+z*z;
    	return (mu_zero/c_light)*(1/(16*M_PI*M_PI))*(1/(r2))*(edpX_dd[i]*edpX_dd[i]+edpY_dd[i]*edpY_dd[i]-(1/(r2*r2))*(x*edpX_dd[i]+y*edpY_dd[i])*(x*edpX_dd[i]+y*edpY_dd[i]));
    };

    auto gr_S = new TGraph();
    int Np = 300;
    double Raio = 10;
    Double_t theta[Np];
   	Double_t radius[Np];

   	for(int j=0; j<Np; j++){
	    theta[j] = 2*M_PI*j/Np;
		radius[j] = Poyting(11,Raio*cos(0)*sin(theta[j]),Raio*sin(0)*sin(theta[j]),Raio*cos(theta[j]));
		gr_S->SetPoint(j,theta[j],Poyting(11,cos(M_PI/4)*sin(theta[j]),sin(M_PI/4)*sin(theta[j]),cos(theta[j])));
	}

   	auto grPolar = new TGraphPolar(Np, theta, radius);
   	grPolar->SetLineColor(2);
    grPolar->SetLineWidth(3);

    c1->Clear();
    grPolar->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta_polar.pdf");

    c1->Clear();
    gr_S->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta.pdf");

   	/*for (int i=0; i<30; i++){
        for(int j=0; j<Np; j++){
	    	theta[j] = 2*M_PI*j/Np;
			radius[j] = Poyting(i,cos(1)*sin(theta[j]),sin(1)*sin(theta[j]),cos(theta[j]));
	    }
	    grPolar = new TGraphPolar(Np, theta, radius);
	    grPolar->SetLineColor(2);
    	grPolar->SetLineWidth(3);

	    c1->Clear();
		grPolar->Draw("AOL");
    	cout << i << endl;
        c1->Print("simulation.gif+5");
    }
	c1->Print("simulation.gif++"); */


	printf("Time taken: %.7fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC); 
	return 0;
}