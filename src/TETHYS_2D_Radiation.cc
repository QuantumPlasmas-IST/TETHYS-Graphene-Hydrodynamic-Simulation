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

#include "../src/Library_Rui/Spline3Interpolator.h"

double c_light = 2.99792458e8;
double mu_zero = 4*M_PI*1e-7;
double epsilon_zero = 8.62607015e-34;

double vF = 10.5;

double D_Conductor = 1;
double Radius_Conductor = 0.;

using namespace std;


bool InConductor(double x, double y){
    if(x*x + y*y < (Radius_Conductor*Radius_Conductor))
        return 1;
    else
        return 0;
}


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
    int Nc = 27;

    vector<double> time;
    vector<double> CurS;
    vector<double> DipX;
    vector<double> DipX_dd;
    vector<double> DipY;
    vector<double> DipY_dd;
    vector<double> QuadXX;
    vector<double> QuadXX_ddd;
    vector<double> QuadXY;
    vector<double> QuadXY_ddd;
    vector<double> QuadYY;
    vector<double> QuadYY_ddd;

    // Q_ij = [Qxx Qxy Qxz] => simetric => [Qxx Qxy Qxz] => no trace=>  [Qxx Qxy    Qxz    ] => 2D (z=0)  => [Qxx  Qxy     0     ] (3 ind. comp.)
    //        [Qyx Qyy Qyz] (9 ind. comp.) [Qxy Qyy Qyz] (6 ind. comp.) [Qxy Qyy    Qyz    ] (5 ind. comp.)  [Qxy  Qyy     0     ]
    //        [Qzx Qzy Qzz]                [Qxy Qyz Qzz]                [Qxz Qyz -(Qxx+Qyy)]                 [ 0    0  -(Qxx+Qyy)]

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
            if(j==6){
                if (stringstream(temp_s) >> temp_d);
                CurS.push_back(temp_d);
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
            if(j==18){
                if (stringstream(temp_s) >> temp_d);
                QuadXX_ddd.push_back(temp_d);
            }
            if(j==19){
                if (stringstream(temp_s) >> temp_d);
                QuadXY.push_back(temp_d);
            }
            if(j==22){
                if (stringstream(temp_s) >> temp_d);
                QuadXY_ddd.push_back(temp_d);
            }
            if(j==23){
                if (stringstream(temp_s) >> temp_d);
                QuadYY.push_back(temp_d);
            }
            if(j==26){
                if (stringstream(temp_s) >> temp_d);
                QuadYY_ddd.push_back(temp_d);
            }
        }
    }

    infile1.close();
    cout << "[Completed reading of " << namefile << " ]" << endl;


    // Interpolation //
    vector<pair<double,double>> DipX_dd_t;
    vector<pair<double,double>> DipY_dd_t;
    vector<pair<double,double>> QuadXX_ddd_t;
    vector<pair<double,double>> QuadXY_ddd_t;
    vector<pair<double,double>> QuadYY_ddd_t;

    pair<double,double> temp_pair;
    for(int i=0; i<Nl; i++){
        temp_pair.first = time[i];

        temp_pair.second = DipX_dd[i];
        DipX_dd_t.push_back(temp_pair);
        temp_pair.second = DipY_dd[i];
        DipY_dd_t.push_back(temp_pair);
        temp_pair.second = QuadXX_ddd[i];
        QuadXX_ddd_t.push_back(temp_pair);
        temp_pair.second = QuadXY_ddd[i];
        QuadXY_ddd_t.push_back(temp_pair);
        temp_pair.second = QuadYY_ddd[i];
        QuadYY_ddd_t.push_back(temp_pair);
    }

    Spline3Interpolator SDipX_dd_t(DipX_dd_t);
    Spline3Interpolator SDipY_dd_t(DipY_dd_t);
    Spline3Interpolator SQuadXX_ddd_t(QuadXX_ddd_t);
    Spline3Interpolator SQuadXY_ddd_t(QuadXY_ddd_t);
    Spline3Interpolator SQuadYY_ddd_t(QuadYY_ddd_t);
    //*/


    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);

    /* // Dipoles and Quadrupole moment graphs //
    auto gr_CurS = new TGraph();

	auto gr_DipX = new TGraph();
	auto gr_DipY = new TGraph();
    auto gr_DipX_dd = new TGraph();
    auto gr_DipY_dd = new TGraph();

    auto gr_QuadXX = new TGraph();
    auto gr_QuadXY = new TGraph();
    auto gr_QuadYY = new TGraph();
    auto gr_QuadXX_ddd = new TGraph();
    auto gr_QuadXY_ddd = new TGraph();
    auto gr_QuadYY_ddd = new TGraph();

    for(int i=0; i<Nl; i++){
        gr_CurS->SetPoint(i,time[i],CurS[i]);

		gr_DipX->SetPoint(i,time[i],DipX[i]);
		gr_DipY->SetPoint(i,time[i],DipY[i]);
        gr_DipX_dd->SetPoint(i,time[i],DipX_dd[i]);
        gr_DipY_dd->SetPoint(i,time[i],DipY_dd[i]);

        gr_QuadXX->SetPoint(i,time[i],QuadXX[i]);
        gr_QuadXY->SetPoint(i,time[i],QuadXY[i]);
        gr_QuadYY->SetPoint(i,time[i],QuadYY[i]);
        //gr_QuadXX_ddd->SetPoint(i,time[i],QuadXX_ddd[i]);
        gr_QuadXY_ddd->SetPoint(i,time[i],QuadXY_ddd[i]);
        gr_QuadYY_ddd->SetPoint(i,time[i],QuadYY_ddd[i]);
    }
    for(int i=0; i<10000; i++)
        gr_QuadXX_ddd->SetPoint(i,(double)2*i/10000, SQuadXX_ddd_t.Interpolate((double)2*i/10000));

	c1->cd();   gr_CurS->Draw();c1->SaveAs("./Files_Images_PIC/CurS_t.pdf");

    c1->Clear();gr_DipX->Draw();c1->SaveAs("./Files_Images_PIC/DipX_t.pdf");
	c1->Clear();gr_DipY->Draw();c1->SaveAs("./Files_Images_PIC/DipY_t.pdf");
    c1->Clear();gr_DipX_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipX_dd_t.pdf");
    c1->Clear();gr_DipY_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipY_dd_t.pdf");

    c1->Clear();gr_QuadXX->Draw();c1->SaveAs("./Files_Images_PIC/QuadXX_t.pdf");
    c1->Clear();gr_QuadXY->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_t.pdf");
    c1->Clear();gr_QuadYY->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_t.pdf");
    c1->Clear();gr_QuadXX_ddd->Draw("AP");c1->SaveAs("./Files_Images_PIC/QuadXX_ddd_t.pdf");
    c1->Clear();gr_QuadXY_ddd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_ddd_t.pdf");
    c1->Clear();gr_QuadYY_ddd->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_ddd_t.pdf");
    //*/


    // Field Functions //
    auto E_Dip = [&](double t, double x, double y, double z){ //X component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double array[] = {(x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*x/r2-SDipX_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*y/r2-SDipY_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*z/r2};
        Vec E_Dip(3,array);
        return (1/(4*M_PI*r))*E_Dip;
    };

    auto H_Dip = [&](double t, double x, double y, double z){ //X component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double array_r[] = {x,y,z};
        Vec r_vec(3,array_r);

        double array_p[] = {SDipX_dd_t.Interpolate(t_r),SDipY_dd_t.Interpolate(t_r),0};
        Vec p_dd(3,array_p);

        return (-1/(4*M_PI*r2))*r_vec.ex(p_dd);
    };

    auto E_Quad = [&](double t, double x, double y, double z){ //X component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

        double array[] = {(x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z};
        Vec E_Quad(3,array);

        return (1/(24*M_PI*300*vF*r))*E_Quad;
    };

    auto H_Quad = [&](double t, double x, double y, double z){ //X component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

        double array_r[] = {x,y,z};
        Vec r_vec(3,array_r);

        double array_Q[] = {Q_vec_X,Q_vec_Y,Q_vec_Z};
        Vec Q_ddd(3,array_Q);

        return (-1/(24*M_PI*300*vF*r2))*r_vec.ex(Q_ddd);
    };

    auto Poynting = [&](double t, double x, double y, double z){ //X component of total Poynting vector // time is in i
        if(InConductor(x*D_Conductor/(z+2*D_Conductor),y*D_Conductor/(z+2*D_Conductor)) == 0){
            return (E_Dip(t,x,y,z)+E_Quad(t,x,y,z)).ex(H_Dip(t,x,y,z)+H_Quad(t,x,y,z));
        }else{
            return (E_Dip(t,x,y,z)+E_Quad(t,x,y,z)-E_Dip(t,x,y,z+2*D_Conductor)-E_Quad(t,x,y,z+2*D_Conductor)).ex(H_Dip(t,x,y,z)+H_Quad(t,x,y,z)-H_Dip(t,x,y,z+2*D_Conductor)-H_Quad(t,x,y,z+2*D_Conductor));
        }
    };
    /////


    int Np = 1000;
    double R = 0.1;
    Double_t theta[Np];
   	Double_t radius[Np];

   	for(int j=0; j<Np; j++){
	    theta[j] = 2*M_PI*j/Np;
		radius[j] = Poynting(0.5,R*sin(theta[j]),0,R*cos(theta[j])).mod();
	}
    auto grPolar = new TGraphPolar(Np, theta, radius);
    grPolar->SetLineColor(2);
    grPolar->SetLineWidth(3);

    c1->Clear();
    grPolar->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta_polar.pdf");

    auto gr3d = new TGraph2D(150*150);

    double theta_3d;
    double phi_3d;
    double intensity;
    double distance = 3050;
    for (int i=0; i<150; i++){
        for(int j=0; j<150; j++){
            theta_3d = 2*M_PI*i/150;
            phi_3d = M_PI*j/150;
            intensity = distance*distance*Poynting(1.8,distance*cos(theta_3d)*sin(phi_3d),distance*sin(theta_3d)*sin(phi_3d),distance*cos(phi_3d)).mod();
            //if(distance*cos(phi_3d) < -D_Conductor) continue;
            gr3d->SetPoint(i*150+j, intensity*cos(theta_3d)*sin(phi_3d), intensity*sin(theta_3d)*sin(phi_3d), intensity*cos(phi_3d));
        }
    }

    c1->Clear();
    gr3d->Draw("PCOL Fi");
    c1->SaveAs("Files_Images_PIC/S_3D.pdf");
    delete gr3d;


    auto grS_t_D = new TGraph2D(200*200);

    distance = 3000;
    double phase;
    for (int i=0; i<200; i++){
        for(int j=0; j<200; j++){
            phase = 1.5+0.25*j/200;
            D_Conductor = 500.*i/200;
            Radius_Conductor = D_Conductor;

            grS_t_D->SetPoint(i*200+j, phase, D_Conductor, distance*distance*Poynting(phase,0,0,distance).mod());
        }
    }
    c1->Clear();
    grS_t_D->SetTitle("; time; Conductor distance");
    grS_t_D->Draw("surf3");
    c1->SaveAs("Files_Images_PIC/S_3D_t_D.pdf");


	cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    cout << "[1A\033[2K\033[1;32mDONE!\033[0m\n";
    cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}

/* //gif3d_t//
    for (int k=0; k<60; k++){
        auto gr3d_gif = new TGraph2D(150*150);

        for (int i=0; i<150; i++){
            for(int j=0; j<150; j++){
                theta_3d = 2*M_PI*i/150;
                phi_3d = M_PI*j/150;
                intensity = distance*distance*Poynting(1.5+0.25*k/60,distance*cos(theta_3d)*sin(phi_3d),distance*sin(theta_3d)*sin(phi_3d),distance*cos(phi_3d)).mod();
                //if(distance*cos(phi_3d) < -D_Conductor) continue;
                gr3d_gif->SetPoint(i*150+j, intensity*cos(theta_3d)*sin(phi_3d), intensity*sin(theta_3d)*sin(phi_3d), intensity*cos(phi_3d));
            }
        }
        c1->Clear();
        gr3d_gif->Draw("PCOL Fi");
        c1->Print("Files_Images_PIC/S_3D_t.gif+15");
        cout << k << endl;
        delete gr3d_gif;
    }
    c1->Print("Files_Images_PIC/S_3D_t.gif++");
    //

/* //gif3d_D_conductor//
for (int k=0; k<60; k++){
        auto gr3d_gif2 = new TGraph2D(150*150);

        D_Conductor = 400.*k/60;
        Radius_Conductor = D_Conductor;
        for (int i=0; i<150; i++){
            for(int j=0; j<150; j++){
                theta_3d = 2*M_PI*i/150;
                phi_3d = M_PI*j/150;
                intensity = distance*distance*Poynting(1.5,distance*cos(theta_3d)*sin(phi_3d),distance*sin(theta_3d)*sin(phi_3d),distance*cos(phi_3d)).mod();
                if(distance*cos(phi_3d) < -D_Conductor) continue;
                gr3d_gif2->SetPoint(i*150+j, intensity*cos(theta_3d)*sin(phi_3d), intensity*sin(theta_3d)*sin(phi_3d), intensity*cos(phi_3d));
            }
        }
        c1->Clear();
        gr3d_gif2->Draw("PCOL Fi");
        c1->Print("Files_Images_PIC/S_3D_D_Conductor.gif+15");
        cout << k << endl;
        delete gr3d_gif2;
    }
    c1->Print("Files_Images_PIC/S_3D_D_Conductor.gif++");
    */


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




    /*/gif - Raio do Condutor a aumentar//
    for (int i=0; i<30; i++){
        D_Conductor = 0.1;
        Radius_Conductor = 0.001*pow(10,(double)i/8);
        for(int j=0; j<Np; j++){
            theta[j] = 2*M_PI*j/Np;
            radius[j] = Poynting_Module(11,R*sin(theta[j]),0,R*cos(theta[j]));
        }
        auto grPolar = new TGraphPolar(Np, theta, radius);
        grPolar->SetLineColor(2);
        grPolar->SetLineWidth(3);
        c1->Clear();
        grPolar->Draw("AOL");
        cout << i << endl;
        c1->Print("Files_Images_PIC/S_Conductor_Radius.gif+15");
    }
    c1->Print("Files_Images_PIC/S_Conductor_Radius.gif++");
    ///*/

    /*/gif - Distancia ao Condutor a diminuir//
    for (int i=0; i<30; i++){
        D_Conductor = pow(0.1,((double)i-10)/3);
        Radius_Conductor = 0.1;
        for(int j=0; j<Np; j++){
            theta[j] = 2*M_PI*j/Np;
            radius[j] = Poynting_Module(11,R*sin(theta[j]),0,R*cos(theta[j]));
        }
        auto grPolar = new TGraphPolar(Np, theta, radius);
        grPolar->SetLineColor(2);
        grPolar->SetLineWidth(3);
        c1->Clear();
        grPolar->Draw("AOL");
        cout << i << endl;
        c1->Print("Files_Images_PIC/S_Conductor_Distance.gif+15");
    }
    c1->Print("Files_Images_PIC/S_Conductor_Distance.gif++");
    ///*/