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

double D_Conductor = 0.4;
double Radius_Conductor = 0.2;

using namespace std;


double Cross_Product_X(double a1, double a2, double a3, double b1, double b2, double b3){
    return a2*b3 - a3*b2;
}
double Cross_Product_Y(double a1, double a2, double a3, double b1, double b2, double b3){
    return a3*b1 - a1*b3;
}
double Cross_Product_Z(double a1, double a2, double a3, double b1, double b2, double b3){
    return a1*b2 - a2*b1;
}

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
    cout << "[Completed reading of " << namefile << " ]" << endl;


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
    c1->Clear();gr_QuadXX_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXX_dd_t.pdf");
    c1->Clear();gr_QuadXY_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_dd_t.pdf");
    c1->Clear();gr_QuadYY_dd->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_dd_t.pdf");

    //Dipole Functions//
    auto E_x_Dip = [&](int i, double x, double y, double z){ //X component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        return (mu_zero/(4*M_PI*r))*((x*DipX_dd[i]+y*DipY_dd[i])*x/r2-DipX_dd[i]);
    };
    auto E_y_Dip = [&](int i, double x, double y, double z){ //Y component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        return (mu_zero/(4*M_PI*r))*((x*DipX_dd[i]+y*DipY_dd[i])*y/r2-DipY_dd[i]);
    };
    auto E_z_Dip = [&](int i, double x, double y, double z){ //Z component of electric field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        return (mu_zero/(4*M_PI*r))*((x*DipX_dd[i]+y*DipY_dd[i])*z/r2);
    };

    auto H_x_Dip = [&](int i, double x, double y, double z){ //X component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        return -Cross_Product_X(x,y,z,DipX_dd[i],DipY_dd[i],0)/(4*M_PI*c_light*r2);
    };
    auto H_y_Dip = [&](int i, double x, double y, double z){ //Y component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        return -Cross_Product_Y(x,y,z,DipX_dd[i],DipY_dd[i],0)/(4*M_PI*c_light*r2);
    };
    auto H_z_Dip = [&](int i, double x, double y, double z){ //Z component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        return -Cross_Product_Z(x,y,z,DipX_dd[i],DipY_dd[i],0)/(4*M_PI*c_light*r2);
    };
    /////

    // Q_ij = [Qxx Qxy Qxz] => simetric => [Qxx Qxy Qxz] => no trace=>  [Qxx Qxy    Qxz    ] => 2D (z=0)  => [Qxx  Qxy     0     ] (3 ind. comp.)
    //        [Qyx Qyy Qyz] (9 ind. comp.) [Qxy Qyy Qyz] (6 ind. comp.) [Qxy Qyy    Qyz    ] (5 ind. comp.)  [Qxy  Qyy     0     ]
    //        [Qzx Qzy Qzz]                [Qxy Qyz Qzz]                [Qxz Qyz -(Qxx+Qyy)]                 [ 0    0  -(Qxx+Qyy)]

    //Quadrupole Functions//
    auto E_x_Quad = [&](int i, double x, double y, double z){ //X component of electric field of Quadrupole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return (mu_zero/(24*M_PI*c_light*r))*((x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X);
    };
    auto E_y_Quad = [&](int i, double x, double y, double z){ //Y component of electric field of Quadrupole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return (mu_zero/(24*M_PI*c_light*r))*((x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y);
    };
    auto E_z_Quad = [&](int i, double x, double y, double z){ //Z component of electric field of Quadrupole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return (mu_zero/(24*M_PI*c_light*r))*((x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z);
    };

    auto H_x_Quad = [&](int i, double x, double y, double z){ //X component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return Cross_Product_X(x,y,z,Q_vec_X,Q_vec_Y,Q_vec_Z)/(24*M_PI*c_light*c_light*r2);
    };
    auto H_y_Quad = [&](int i, double x, double y, double z){ //Y component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return Cross_Product_Y(x,y,z,Q_vec_X,Q_vec_Y,Q_vec_Z)/(24*M_PI*c_light*c_light*r2);
    };
    auto H_z_Quad = [&](int i, double x, double y, double z){ //Z component of H field of Dipole // time is in i
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double Q_vec_X = QuadXX_dd[i]*x + QuadXY_dd[i]*y;
        double Q_vec_Y = QuadXY_dd[i]*x + QuadYY_dd[i]*y;
        double Q_vec_Z = -(QuadXX_dd[i]+QuadYY_dd[i])*z;
        return Cross_Product_Z(x,y,z,Q_vec_X,Q_vec_Y,Q_vec_Z)/(24*M_PI*c_light*c_light*r2);
    };
    /////

    //Poynting//
    auto Poynting_X = [&](int i, double x, double y, double z){ //X component of total Poynting vector // time is in i
        if(InConductor(x*D_Conductor/(z+2*D_Conductor),y*D_Conductor/(z+2*D_Conductor)) == 0){
            return Cross_Product_X(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z),E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z),E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z),H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z),H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z));
        }else{
            return Cross_Product_X(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z)-E_x_Dip(i,x,y,z)-E_x_Quad(i,x,y,z)
                                  ,E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z)-E_y_Dip(i,x,y,z)-E_y_Quad(i,x,y,z)
                                  ,E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)-E_z_Dip(i,x,y,z)-E_z_Quad(i,x,y,z)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z)-H_x_Dip(i,x,y,z)-H_x_Quad(i,x,y,z)
                                  ,H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z)-H_y_Dip(i,x,y,z)-H_y_Quad(i,x,y,z)
                                  ,H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z)-H_z_Dip(i,x,y,z)-H_z_Quad(i,x,y,z));
        }
    };
    auto Poynting_Y = [&](int i, double x, double y, double z){ //X component of total Poynting vector // time is in i
        if(InConductor(x*D_Conductor/(z+2*D_Conductor),y*D_Conductor/(z+2*D_Conductor)) == 0){
            return Cross_Product_Y(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z),E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z),E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z),H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z),H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z));
        }else{
            return Cross_Product_Y(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z)-E_x_Dip(i,x,y,z)-E_x_Quad(i,x,y,z)
                                  ,E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z)-E_y_Dip(i,x,y,z)-E_y_Quad(i,x,y,z)
                                  ,E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)-E_z_Dip(i,x,y,z)-E_z_Quad(i,x,y,z)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z)-H_x_Dip(i,x,y,z)-H_x_Quad(i,x,y,z)
                                  ,H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z)-H_y_Dip(i,x,y,z)-H_y_Quad(i,x,y,z)
                                  ,H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z)-H_z_Dip(i,x,y,z)-H_z_Quad(i,x,y,z));
        }
    };
    auto Poynting_Z = [&](int i, double x, double y, double z){ //X component of total Poynting vector // time is in i
        if(InConductor(x*D_Conductor/(z+2*D_Conductor),y*D_Conductor/(z+2*D_Conductor)) == 0){
            return Cross_Product_Z(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z),E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z),E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z),H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z),H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z));
        }else{
            return Cross_Product_Z(E_x_Dip(i,x,y,z)+E_x_Quad(i,x,y,z)-E_x_Dip(i,x,y,z+2*D_Conductor)-E_x_Quad(i,x,y,z+2*D_Conductor)
                                  ,E_y_Dip(i,x,y,z)+E_y_Quad(i,x,y,z)-E_y_Dip(i,x,y,z+2*D_Conductor)-E_y_Quad(i,x,y,z+2*D_Conductor)
                                  ,E_z_Dip(i,x,y,z)+E_z_Quad(i,x,y,z)-E_z_Dip(i,x,y,z+2*D_Conductor)-E_z_Quad(i,x,y,z+2*D_Conductor)
                                  ,H_x_Dip(i,x,y,z)+H_x_Quad(i,x,y,z)-H_x_Dip(i,x,y,z+2*D_Conductor)-H_x_Quad(i,x,y,z+2*D_Conductor)
                                  ,H_y_Dip(i,x,y,z)+H_y_Quad(i,x,y,z)-H_y_Dip(i,x,y,z+2*D_Conductor)-H_y_Quad(i,x,y,z+2*D_Conductor)
                                  ,H_z_Dip(i,x,y,z)+H_z_Quad(i,x,y,z)-H_z_Dip(i,x,y,z+2*D_Conductor)-H_z_Quad(i,x,y,z+2*D_Conductor));
        }
    };
    auto Poynting_Module = [&](int i, double x, double y, double z){ //X component of total Poynting vector // time is in i
        return sqrt(Poynting_X(i,x,y,z)*Poynting_X(i,x,y,z)+Poynting_Y(i,x,y,z)*Poynting_Y(i,x,y,z)+Poynting_Z(i,x,y,z)*Poynting_Z(i,x,y,z));
    };
    ////


    int Np = 1000;
    double R = 0.1;
    Double_t theta[Np];
   	Double_t radius[Np];

   	for(int j=0; j<Np; j++){
	    theta[j] = 2*M_PI*j/Np;
		radius[j] = Poynting_Module(11,R*sin(theta[j]),0,R*cos(theta[j]));
	}
    auto grPolar = new TGraphPolar(Np, theta, radius);
    grPolar->SetLineColor(2);
    grPolar->SetLineWidth(3);

    c1->Clear();
    grPolar->Draw("AOL");
    c1->SaveAs("Files_Images_PIC/S_theta_polar.pdf");

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


    auto gr3d = new TGraph2D(150*150);
    gr3d->SetName("gr1_name");

    double theta_3d;
    double phi_3d;
    double intensity;
    for (int i=0; i<150; i++){
        for(int j=0; j<150; j++){
            theta_3d = 2*M_PI*i/150;
            phi_3d = M_PI*j/150;
            intensity = 1e14*Poynting_Module(11,cos(theta_3d)*sin(phi_3d),sin(theta_3d)*sin(phi_3d),cos(phi_3d));
            if(cos(phi_3d) < -D_Conductor) continue;
            gr3d->SetPoint(i*150+j, intensity*cos(theta_3d)*sin(phi_3d), intensity*sin(theta_3d)*sin(phi_3d), intensity*cos(phi_3d));
        }
    }

    c1->Clear();
    gr3d->SetTitle("Radiation_Mirrored_Dipole+Quadrupole");
    gr3d->SetMarkerColor(kBlue);
    gr3d->Draw("PCOL Fi");
    c1->SaveAs("Files_Images_PIC/S_3D.pdf");

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