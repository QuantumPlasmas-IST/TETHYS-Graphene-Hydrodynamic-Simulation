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

double epsilon_r = 3;
double mu_r = 1.5;

double vF = 10.5;

double n_index = sqrt(epsilon_r*mu_r);
double D_Conductor = 3;
double Radius_Conductor = 3;
double D_dieletric = 0.1;
double distance_measure = 3000;
double d_xy, global_z;

using namespace std;

double NewtonsMethod(double(*func)(double), double x1=1.56, double x2=1.57, double error_min=1e-5, int iter_max=40){
    int i = 0;
    double fx1 = 0., fx2 = 0., x3 = 0.;
    while((fabs(x1-x2) >= error_min || fabs(fx2) > error_min) && i < iter_max){
        fx2 = func(x2);
        fx1 = func(x1);
        if(fx2-fx1 != 0) x3 = x2 - (fx2*(x2-x1)/(fx2-fx1));
        else{
            cout << "We have found a stationary point at x = " << x2 << ". Ending Newton's method here." << endl;
            return x2;
        }
        if(x3>1.6) x3=M_PI-1e-3;
        if(x3<0) x3=M_PI-2e-3;
        x1 = x2;
        x2 = x3;
        i += 1;
    }
    if(i == iter_max) cout << "[NewtonsMethod] Iteration limit reached" << endl;
    return x2;
}

double Theta_refraction_equation(double theta_r){
    return (global_z+2*D_dieletric)*tan(theta_r) + 2*(D_Conductor-D_dieletric)*tan(asin(sin(theta_r)/n_index)) - d_xy;
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
    int Nc = 27;

    vector<double> time, CurS;
    vector<double> DipX, DipX_dd, DipY, DipY_dd;
    vector<double> QuadXX, QuadXX_ddd, QuadXY, QuadXY_ddd, QuadYY, QuadYY_ddd;

    int Nl = 0; //number of lines in file
    while(!infile1.eof()){
        Nl++;
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
    Nl=Nl-2;

    infile1.close();
    cout << "[Completed reading of " << namefile << " ]" << "\n" << endl;

    int Nl_0 = 1700;
    cout << "The last time point is " << time[Nl] << " at line " << Nl << endl;
    cout << "We will work in the following time interval [" << time[Nl_0] << " , " << time[Nl] << "]\n"<< endl;


    // Interpolation //
        vector<pair<double,double>> DipX_dd_t;
        vector<pair<double,double>> DipY_dd_t;
        vector<pair<double,double>> QuadXX_ddd_t;
        vector<pair<double,double>> QuadXY_ddd_t;
        vector<pair<double,double>> QuadYY_ddd_t;


        pair<double,double> temp_pair;
        for(int i=Nl_0; i<Nl; i++){
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

        cout << "Cubic Spline Interpolation Complete" << endl;
        //*/


    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);

    // Dipoles and Quadrupole moment graphs //

    /*/ Original points graphs
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

        for(int i=Nl_0; i<Nl; i++){
            gr_CurS->SetPoint(i-Nl_0,time[i],CurS[i]);

    		gr_DipX->SetPoint(i-Nl_0,time[i],DipX[i]);
    		gr_DipY->SetPoint(i-Nl_0,time[i],DipY[i]);
            gr_DipX_dd->SetPoint(i-Nl_0,time[i],DipX_dd[i]);
            gr_DipY_dd->SetPoint(i-Nl_0,time[i],DipY_dd[i]);

            gr_QuadXX->SetPoint(i-Nl_0,time[i],QuadXX[i]);
            gr_QuadXY->SetPoint(i-Nl_0,time[i],QuadXY[i]);
            gr_QuadYY->SetPoint(i-Nl_0,time[i],QuadYY[i]);
            gr_QuadXX_ddd->SetPoint(i-Nl_0,time[i],QuadXX_ddd[i]);
            gr_QuadXY_ddd->SetPoint(i-Nl_0,time[i],QuadXY_ddd[i]);
            gr_QuadYY_ddd->SetPoint(i-Nl_0,time[i],QuadYY_ddd[i]);
        }
        c1->cd();   gr_CurS->Draw();c1->SaveAs("./Files_Images_PIC/CurS_t.pdf");

        c1->Clear(); gr_DipX->Draw();c1->SaveAs("./Files_Images_PIC/DipX_t.pdf");
        c1->Clear(); gr_DipY->Draw();c1->SaveAs("./Files_Images_PIC/DipY_t.pdf");
        c1->Clear(); gr_DipX_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipX_dd_t.pdf");
        c1->Clear(); gr_DipY_dd->Draw();c1->SaveAs("./Files_Images_PIC/DipY_dd_t.pdf");

        c1->Clear(); gr_QuadXX->Draw();c1->SaveAs("./Files_Images_PIC/QuadXX_t.pdf");
        c1->Clear(); gr_QuadXY->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_t.pdf");
        c1->Clear(); gr_QuadYY->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_t.pdf");
        c1->Clear(); gr_QuadXX_ddd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXX_ddd_t.pdf");
        c1->Clear(); gr_QuadXY_ddd->Draw();c1->SaveAs("./Files_Images_PIC/QuadXY_ddd_t.pdf");
        c1->Clear(); gr_QuadYY_ddd->Draw();c1->SaveAs("./Files_Images_PIC/QuadYY_ddd_t.pdf");
        //*/

    /*/ Interpolated points graphs
        auto grI_CurS = new TGraph();

        auto grI_DipX_dd = new TGraph();
        auto grI_DipY_dd = new TGraph();

        auto grI_QuadXX_ddd = new TGraph();
        auto grI_QuadXY_ddd = new TGraph();
        auto grI_QuadYY_ddd = new TGraph();

        int Np_interpolation = 10000;
        double t_interpolation = time[Nl_0];
        for(int i=0; i<Np_interpolation; i++){
            grI_DipX_dd->SetPoint(i,t_interpolation,SDipX_dd_t.Interpolate(t_interpolation));
            grI_DipY_dd->SetPoint(i,t_interpolation,SDipY_dd_t.Interpolate(t_interpolation));

            grI_QuadXX_ddd->SetPoint(i,t_interpolation,SQuadXX_ddd_t.Interpolate(t_interpolation));
            grI_QuadXY_ddd->SetPoint(i,t_interpolation,SQuadXY_ddd_t.Interpolate(t_interpolation));
            grI_QuadYY_ddd->SetPoint(i,t_interpolation,SQuadYY_ddd_t.Interpolate(t_interpolation));

            t_interpolation += (time[Nl]-time[Nl_0])/Np_interpolation;
        }
        c1->Clear(); grI_DipX_dd->Draw("AP");c1->SaveAs("./Files_Images_PIC/IDipX_dd_t.pdf");
        c1->Clear(); grI_DipY_dd->Draw("AP");c1->SaveAs("./Files_Images_PIC/IDipY_dd_t.pdf");

        c1->Clear(); grI_QuadXX_ddd->Draw("AP");c1->SaveAs("./Files_Images_PIC/IQuadXX_ddd_t.pdf");
        c1->Clear(); grI_QuadXY_ddd->Draw("AP");c1->SaveAs("./Files_Images_PIC/IQuadXY_ddd_t.pdf");
        c1->Clear(); grI_QuadYY_ddd->Draw("AP");c1->SaveAs("./Files_Images_PIC/IQuadYY_ddd_t.pdf");
        //*/


    // Field Functions //
    auto E_field = [&](double t, double x, double y, double z, double theta, double theta_i){
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double array[] = {(x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*x/r2-SDipX_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*y/r2-SDipY_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*z/r2};
        Vec E_Dip(3,array);

        double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

        double array2[] = {(x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z};
        Vec E_Quad(3,array2);

        Vec E_field = (1/(4*M_PI*r))*E_Dip + 0*(1/(24*M_PI*300*vF*r))*E_Quad;

        double theta_r = asin(sin(theta_i)/n_index);

        double r_T, r_ll;
        if(theta_i<1e-3){
            r_T = (1-n_index)/(1+n_index);
            r_ll = r_T;
        }else{
            r_T  = -sin(theta_i-theta_r)/sin(theta_i+theta_r); //perpendicular reflected coeficient
            r_ll = tan(theta_i-theta_r)/tan(theta_i+theta_r);  //parallel refrected coeficient
        }

        double array_T[] = {E_field[0]*sin(theta), E_field[1]*cos(theta), 0};
        Vec E_field_T(3,array_T);

        return E_field+(E_field_T*r_T + (E_field-E_field_T)*r_ll);
    };

    auto E_field_Image = [&](double t, double x, double y, double z, double theta, double theta_r){
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);

        double theta_i = asin(sin(theta_r)/n_index);

        double r_vacuum;
        double r_dieletric = 2*(D_Conductor - D_dieletric)/cos(theta_i);
        if(theta_r>1.5707) d_xy - sin(theta_i)*r_dieletric;
        else r_vacuum = (z + 2*D_dieletric)/cos(theta_r);
        double t_r = t - (n_index*r_dieletric+r_vacuum)/(300*vF);

        double array[] = {(x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*x/r2-SDipX_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*y/r2-SDipY_dd_t.Interpolate(t_r),
                          (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*z/r2};
        Vec E_Dip(3,array);

        double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
        double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

        double array2[] = {(x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y,
                          (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z};
        Vec E_Quad(3,array2);

        Vec E_field = (1/(4*M_PI*r))*E_Dip + (1/(24*M_PI*300*vF*r))*E_Quad;

        double t_T, t_ll;
        if(theta_i<1e-3){
            t_T = 4*n_index/((n_index+1)*(n_index+1));
            t_ll = t_T;
        }else{
            t_T  = 4*sin(theta_r)*sin(theta_i)*cos(theta_i)*cos(theta_r)/(sin(theta_i+theta_r)*sin(theta_i+theta_r)); //perpendicular transmited coeficient
            t_ll = t_T/(cos(theta_i - theta_r)*cos(theta_i - theta_r));                       //parallel transmited coeficient
        }

        double array_T[] = {E_field[0]*sin(theta), E_field[1]*cos(theta), 0};
        Vec E_field_T(3,array_T); 

        return -(E_field_T*t_T + (E_field-E_field_T)*t_ll);
    };

    auto H_field = [&](double t, double x, double y, double z, double theta, double theta_i){
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);
        double t_r = t-r/(300*vF);

        double array_r[] = {x,y,z};
        Vec r_vec(3,array_r);

        double array_p[] = {SDipX_dd_t.Interpolate(t_r),SDipY_dd_t.Interpolate(t_r),0};
        Vec p_dd(3,array_p);

        double array_Q[] = {(SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r
                            ,(SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r
                            ,-(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r};
        Vec Q_ddd(3,array_Q);

        Vec H_field = (-1/(4*M_PI*r2))*r_vec.ex(p_dd) + 0*(-1/(24*M_PI*300*vF*r2))*r_vec.ex(Q_ddd);

        double theta_r = asin(sin(theta_i)/n_index);

        double r_T, r_ll;
        if(theta_i<1e-3){
            r_T = (1-n_index)/(1+n_index);
            r_ll = r_T;
        }else{
            r_T  = -sin(theta_i-theta_r)/sin(theta_i+theta_r); //perpendicular reflected coeficient
            r_ll = tan(theta_i-theta_r)/tan(theta_i+theta_r);  //parallel refrected coeficient
        }

        double array_T[] = {H_field[0]*sin(theta), H_field[1]*cos(theta), 0};
        Vec H_field_T(3,array_T);

        return H_field+(H_field_T*r_T + (H_field-H_field_T)*r_ll);
    };

    auto H_field_Image = [&](double t, double x, double y, double z, double theta, double theta_r){
        double r2 = x*x+y*y+z*z;
        double r = sqrt(r2);

        double theta_i = asin(sin(theta_r)/n_index);

        double r_vacuum;
        double r_dieletric = 2*(D_Conductor - D_dieletric)/cos(theta_i);
        if(theta_r>1.5707) d_xy - sin(theta_i)*r_dieletric;
        else r_vacuum = (z + 2*D_dieletric)/cos(theta_r);
        double t_r = t - (n_index*r_dieletric+r_vacuum)/(300*vF);

        double array_r[] = {x,y,z};
        Vec r_vec(3,array_r);

        double array_p[] = {SDipX_dd_t.Interpolate(t_r),SDipY_dd_t.Interpolate(t_r),0};
        Vec p_dd(3,array_p);

        double array_Q[] = {(SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r
                            ,(SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r
                            ,-(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r};
        Vec Q_ddd(3,array_Q);

        Vec H_field = (-1/(4*M_PI*r2))*r_vec.ex(p_dd) + (-1/(24*M_PI*300*vF*r2))*r_vec.ex(Q_ddd);

        double t_T, t_ll;
        if(theta_i<1e-3){
            t_T = 4*n_index/((n_index+1)*(n_index+1));
            t_ll = t_T;
        }else{
            t_T  = 4*sin(theta_r)*sin(theta_i)*cos(theta_i)*cos(theta_r)/(sin(theta_i+theta_r)*sin(theta_i+theta_r)); //perpendicular transmited coeficient
            t_ll = t_T/(cos(theta_i - theta_r)*cos(theta_i - theta_r));                       //parallel transmited coeficient
        }

        double array_T[] = {H_field[0]*sin(theta), H_field[1]*cos(theta), 0};
        Vec H_field_T(3,array_T); 

        return -(H_field_T*t_T + (H_field-H_field_T)*t_ll);
    };

    auto Emitter = [&](double t, double O_x, double O_y, double x, double y, double z){
        d_xy = sqrt((x-O_x)*(x-O_x)+(y-O_y)*(y-O_y));
        global_z = z;
        if(n_index < 1+1e-5) n_index = 1+1e-5;

        double theta_r;
        double d_star;
        if(fabs(z) < 1e-3){
            theta_r = M_PI/2;
            d_star = D_Conductor*sin(asin(1/n_index));
        }
        else{
            theta_r = NewtonsMethod(&Theta_refraction_equation);
            d_star = d_xy - (z + D_dieletric)*tan(theta_r) - (D_Conductor - D_dieletric)*tan(asin(sin(theta_r)/n_index));
        }

        double theta;
        if(fabs(x) < 1e-3)
            theta = M_PI/2;
        else 
            theta = fabs(atan((y - O_y)/(x - O_x)));

        double theta_i_original = atan(sqrt(d_xy/z));

        if(InConductor(O_x + d_star*(x-O_x)/d_xy, O_y + d_star*(y-O_y)/d_xy) == 0 && d_xy > 1e-3){
            return make_pair(E_field(t,x-O_x,y-O_y,z,theta,theta_i_original),H_field(t,x-O_x,y-O_y,z,theta,theta_i_original));
        }else{
            return make_pair(E_field(t,x-O_x,y-O_y,z,theta,theta_i_original)+E_field_Image(t,x-O_x,y-O_y,z,theta,theta_r),H_field(t,x-O_x,y-O_y,z,theta,theta_i_original)+H_field_Image(t,x-O_x,y-O_y,z,theta,theta_r));
        }
    };

    auto Poyting = [&](Vec E, Vec H){
        return E.ex(H);
    }; 
    /////

    // 3Dgraph
        auto gr3d = new TGraph2D(500*500);

        double theta_3d;
        double phi_3d;
        double intensity;
        for (int i=0; i<500; i++){
            for(int j=0; j<500; j++){
                theta_3d = 2*M_PI*i/500;
                phi_3d = M_PI*j/500;
                if(cos(phi_3d) < 0) continue;
                auto fields = Emitter(9, 0, 0, distance_measure*cos(theta_3d)*sin(phi_3d),distance_measure*sin(theta_3d)*sin(phi_3d),distance_measure*cos(phi_3d));
                intensity = distance_measure*distance_measure*Poyting(fields.first, fields.second).mod();
                gr3d->SetPoint(i*500+j, intensity*cos(theta_3d)*sin(phi_3d), intensity*sin(theta_3d)*sin(phi_3d), intensity*cos(phi_3d));
            }
        }

        c1->Clear();
        gr3d->Draw("PCOL Fi");
        c1->SaveAs("Files_Images_PIC/S_3D.pdf");
        delete gr3d;
        //*/

    /*/S Integrated over time = f(D_Conductor)
        double delta_t = 2;
        int N_points = 200;
        auto grS_D_integrated = new TGraph(N_points);

        double phase;
        double integral;
        for (int i=0; i<1000; i++){
            D_Conductor = 400.*i/1000;
            Radius_Conductor = D_Conductor + 0.2;

            integral = 0;
            for(int j=0; j<200; j++){
                phase = 9+(double)j*delta_t/N_points;
                auto fields = Emitter(phase,0,0,0,0,distance_measure);
                if(j == 0 || j == N_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fields.first, fields.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fields.first, fields.second).mod();
            }
            grS_D_integrated->SetPoint(i, D_Conductor, integral*delta_t/N_points);
        }
        c1->Clear();
        grS_D_integrated->SetTitle("; Conductor distance_measure; Integrated |S|");
        grS_D_integrated->Draw("AP");
        c1->SaveAs("Files_Images_PIC/S_D_integrated.pdf");
        //*/


    cout << "[1A\033[2K\033[1;32mDONE!\033[0m\n";
    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}






    /*/ Polar graph
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
        //*/

    /* //gif3d_t//
        for (int k=0; k<60; k++){
            auto gr3d_gif = new TGraph2D(150*150);

            for (int i=0; i<150; i++){
                for(int j=0; j<150; j++){
                    theta_3d = 2*M_PI*i/150;
                    phi_3d = M_PI*j/150;
                    intensity = distance_measure*distance_measure*Poynting(1.5+0.25*k/60,distance_measure*cos(theta_3d)*sin(phi_3d),distance_measure*sin(theta_3d)*sin(phi_3d),distance_measure*cos(phi_3d)).mod();
                    //if(distance_measure*cos(phi_3d) < -D_Conductor) continue;
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
                    intensity = distance_measure*distance_measure*Poynting(1.5,distance_measure*cos(theta_3d)*sin(phi_3d),distance_measure*sin(theta_3d)*sin(phi_3d),distance_measure*cos(phi_3d)).mod();
                    if(distance_measure*cos(phi_3d) < -D_Conductor) continue;
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
            c1->Print("Files_Images_PIC/S_Conductor_distance_measure.gif+15");
        }
        c1->Print("Files_Images_PIC/S_Conductor_distance_measure.gif++");
        ///*/