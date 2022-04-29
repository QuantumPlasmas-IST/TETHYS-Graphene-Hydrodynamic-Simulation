#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

//O ROOT só está a fazer gráficos e é temporario para ir vendo os resultados mais rapidamente
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

int main(int argc, char **argv){
    string input_file_name;
    if(argc==2){
        input_file_name = argv[1];
    }else{
        cout <<"Input whole name of electro_2D_[...].dat file to read:\t";
        cin >> input_file_name;
    }
	clock_t tStart = clock();

	string namefile = "build/" + input_file_name;
	ifstream infile1(namefile, ios::in);

    string line;
    double temp_d;
    int Nl = 474; //number of lines in file
    int Nc = 24;

    vector<double> time;
    vector<double> DipX;
    vector<double> DipX_dd;
    vector<double> DipY;
    vector<double> DipY_dd;
    vector<double> QuadXX;
    vector<double> QuadXX_dd;
    vector<double> QuadXY;
    vector<double> QuadXY_dd;
    vector<double> QuadYY;
    vector<double> QuadYY_dd;

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
            	DipX.push_back(temp_d);
            }
            if(j==11){
            	if (stringstream(temp_s) >> temp_d);
            	DipX_dd.push_back(temp_d);
            }
            if(j==12){
            	if (stringstream(temp_s) >> temp_d);
            	DipY.push_back(temp_d);
            }
            if(j==14){
            	if (stringstream(temp_s) >> temp_d);
            	DipY_dd.push_back(temp_d);
            }
            if(j==15){
                if (stringstream(temp_s) >> temp_d);
                QuadXX.push_back(temp_d);
            }
            if(j==17){
                if (stringstream(temp_s) >> temp_d);
                QuadXX_dd.push_back(temp_d);
            }
            if(j==15){
                if (stringstream(temp_s) >> temp_d);
                QuadXX.push_back(temp_d);
            }
            if(j==17){
                if (stringstream(temp_s) >> temp_d);
                QuadXX_dd.push_back(temp_d);
            }
            if(j==18){
                if (stringstream(temp_s) >> temp_d);
                QuadXY.push_back(temp_d);
            }
            if(j==20){
                if (stringstream(temp_s) >> temp_d);
                QuadXY_dd.push_back(temp_d);
            }
            if(j==21){
                if (stringstream(temp_s) >> temp_d);
                QuadYY.push_back(temp_d);
            }
            if(j==23){
                if (stringstream(temp_s) >> temp_d);
                QuadYY_dd.push_back(temp_d);
            }
        }
    }

    infile1.close();
    cout << "[Reading Complete " << namefile << " ]" << endl;


    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
	auto gr_DipX = new TGraph();
	auto gr_DipY = new TGraph();
    auto gr_DipX_dd = new TGraph();
    auto gr_DipY_dd = new TGraph();

    auto gr_QuadXX = new TGraph();
    auto gr_QuadXY = new TGraph();
    auto gr_QuadYY = new TGraph();
    auto gr_QuadXX_dd = new TGraph();
    auto gr_QuadXY_dd = new TGraph();
    auto gr_QuadYY_dd = new TGraph();

    for(int i=0; i<Nl; i++){
		gr_DipX->SetPoint(i,time[i],DipX[i]);
		gr_DipY->SetPoint(i,time[i],DipY[i]);
        gr_DipX_dd->SetPoint(i,time[i],DipX_dd[i]);
        gr_DipY_dd->SetPoint(i,time[i],DipY_dd[i]);

        gr_QuadXX->SetPoint(i,time[i],QuadXX[i]);
        gr_QuadXY->SetPoint(i,time[i],QuadXY[i]);
        gr_QuadYY->SetPoint(i,time[i],QuadYY[i]);
        gr_QuadXX_dd->SetPoint(i,time[i],QuadXX_dd[i]);
        gr_QuadXY_dd->SetPoint(i,time[i],QuadXY_dd[i]);
        gr_QuadYY_dd->SetPoint(i,time[i],QuadYY_dd[i]);
    }

	c1->cd();   gr_DipX->Draw();c1->SaveAs("./Files_Images_PIC/DipX_t.pdf");
	c1->Clear();gr_DipY->Draw();c1->SaveAs("./Files_Images_PIC/DipY_t.pdf");
    c1->Clear();gr_DipX_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipX_dd_t.pdf");
    c1->Clear();gr_DipY_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipY_dd_t.pdf");

    c1->Clear();gr_QuadXX->Draw();c1->SaveAs("./Files_Images_PIC/QuadXX_t.pdf");
    c1->Clear();gr_QuadXY->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_t.pdf");
    c1->Clear();gr_QuadYY->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_t.pdf");
    c1->Clear();gr_QuadXX_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_dd_t.pdf");
    c1->Clear();gr_QuadXY_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_dd_t.pdf");
    c1->Clear();gr_QuadYY_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_dd_t.pdf");


	auto Poyting = [&](int i, double x, double y, double z){ //modulo do vetor Poyting // tempo está no i
		double r2 = x*x+y*y+z*z;
    	return (mu_zero/c_light)*(1/(16*M_PI*M_PI))*(1/(r2))*(DipX_dd[i]*DipX_dd[i]+DipY_dd[i]*DipY_dd[i]-(1/(r2*r2))*(x*DipX_dd[i]+y*DipY_dd[i])*(x*DipX_dd[i]+y*DipY_dd[i]));
    };

    auto gr_S = new TGraph();
    double Raio = 10;
    Double_t theta[Nl];
   	Double_t radius[Nl];

   	for(int j=0; j<Nl; j++){
	    theta[j] = 2*M_PI*j/Nl;
		radius[j] = Poyting(11,Raio*cos(0)*sin(theta[j]),Raio*sin(0)*sin(theta[j]),Raio*cos(theta[j]));
		gr_S->SetPoint(j,theta[j],Poyting(11,cos(M_PI/4)*sin(theta[j]),sin(M_PI/4)*sin(theta[j]),cos(theta[j])));
	}

   	auto grPolar = new TGraphPolar(Nl, theta, radius);
   	grPolar->SetLineColor(2);
    grPolar->SetLineWidth(3);

    c1->Clear();
    grPolar->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta_polar.pdf");

    c1->Clear();
    gr_S->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta.pdf");


	cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    cout << "[1A\033[2K\033[1;32mDONE!\033[0m\n";
    cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}







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