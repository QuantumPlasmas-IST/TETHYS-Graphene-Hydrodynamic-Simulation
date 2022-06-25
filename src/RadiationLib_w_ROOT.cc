/************************************************************************************************\
* 2022 Rui Martins, Pedro Cosme                                                                 *
*                                                                                               *
*                                                                                               *
\************************************************************************************************/

//Version with ROOT
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
//*/

#include "includes/RadiationLib_w_ROOT.h"

#ifndef MAT_PI
#    define MAT_PI 3.14159265358979323846
#endif

#ifndef MAT_EULER
#    define MAT_EULER 2.71828182845905
#endif

using namespace std;

//Constructor
Radiation::Radiation(double n_index_ref, double epsilon_r_input, double Dist_Conductor, double R_Conductor, double Distance_measure){
    if (n_index_ref < (1-1e-3)) cout << "Refractive index must be bigger than 1\n";
    if (epsilon_r_input < (1-1e-3)) cout << "Relative permittivity must be bigger than 1\n";
    if (Dist_Conductor<0) cout << "Distance to conductor must be bigger than 0 \n";
    if (R_Conductor<0) cout << "Radius must not be smaller than 0 \n";
    if (Distance_measure<0) cout << "Measurement distance must be bigger than 0 \n";

    n_index = n_index_ref;
    epsilon_r = epsilon_r_input;
    D_Conductor = Dist_Conductor;
    Radius_Conductor= R_Conductor;
    distance_measure = Distance_measure;

    mu_r = n_index*n_index/epsilon_r;

    if(mu_r < 1) cout << "relative permeability can't be smaller than 1 \n";

    rotation180 = 1;
}


//Getting points
int Radiation::ReadFromElectroFile(const string input_file_name){
    vF = stod(input_file_name.substr(input_file_name.find("vF=")+3,input_file_name.find("vis=")-input_file_name.find("vF=")-3));

    string namefile = input_file_name;
    ifstream infile1(namefile, ios::in);

    string line;
    double temp_d;
    int Nc = 30;

    int Nl = 0; //number of lines in file
    while(!infile1.eof()){
        Nl++;
        getline(infile1, line);
        stringstream ss;
        ss << line;

        for(int j=0; j<Nc; j++){
            if (ss >> temp_d){
                if(j==0) time.push_back(temp_d);
                if(j==6) CurS.push_back(temp_d);
                if(j==9) DipX.push_back(temp_d);
                if(j==11)DipX_dd.push_back(temp_d);
                if(j==12)DipY.push_back(temp_d);
                if(j==14)DipY_dd.push_back(temp_d);
                if(j==15)QuadXX.push_back(temp_d);
                if(j==18)QuadXX_ddd.push_back(temp_d);
                if(j==19)QuadXY.push_back(temp_d);
                if(j==22)QuadXY_ddd.push_back(temp_d);
                if(j==23)QuadYY.push_back(temp_d);
                if(j==26)QuadYY_ddd.push_back(temp_d);
                if(j==27) DipMagZ.push_back(temp_d);
                if(j==29)DipMagZ_dd.push_back(temp_d);
            }
        }
    }

    infile1.close();
    cout << "The last time point is " << time[Nl-2] << " at line " << Nl-2 << endl;
    cout << "[Completed reading of " << namefile << " ]" << "\n" << endl;

    Nl=Nl-2;
    return Nl;
}

void Radiation::Interpolation(int Nl_0, int Nl){
    cout << "We will work in the following time interval [" << time[Nl_0] << " , " << time[Nl] << "]\n"<< endl;

    vector<pair<double,double>> DipX_dd_t;
    vector<pair<double,double>> DipY_dd_t;
    vector<pair<double,double>> QuadXX_ddd_t;
    vector<pair<double,double>> QuadXY_ddd_t;
    vector<pair<double,double>> QuadYY_ddd_t;
    vector<pair<double,double>> DipMagZ_dd_t;

    t_min = time[Nl_0];
    t_max = time[Nl];

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
        temp_pair.second = DipMagZ_dd[i];
        DipMagZ_dd_t.push_back(temp_pair);
    }

    Spline3Interpolator SDipX_dd_t_Dummy(DipX_dd_t);
    Spline3Interpolator SDipY_dd_t_Dummy(DipY_dd_t);
    Spline3Interpolator SQuadXX_ddd_t_Dummy(QuadXX_ddd_t);
    Spline3Interpolator SQuadXY_ddd_t_Dummy(QuadXY_ddd_t);
    Spline3Interpolator SQuadYY_ddd_t_Dummy(QuadYY_ddd_t);
    Spline3Interpolator SDipMagZ_dd_t_Dummy(DipMagZ_dd_t);

    SDipX_dd_t = SDipX_dd_t_Dummy;
    SDipY_dd_t = SDipY_dd_t_Dummy;
    SQuadXX_ddd_t = SQuadXX_ddd_t_Dummy;
    SQuadXY_ddd_t = SQuadXY_ddd_t_Dummy;
    SQuadYY_ddd_t = SQuadYY_ddd_t_Dummy;
    SDipMagZ_dd_t = SDipMagZ_dd_t_Dummy;

    cout << "Cubic Spline Interpolation Complete\n\n";

    int contador = 0;
    double sum_c = 0;
    for (int i=Nl_0; i<Nl; i++){
        sum_c += atan(DipY[i]/DipX[i]);
        contador++;
        if(fabs(DipX[i])<1e-12){
            E_plane_angle = M_PI/2;
            break;
        }
    }
    E_plane_angle = sum_c/contador+M_PI/2;
}


//Printing graphs
void Radiation::PrintGraphElectro(int Nl_0, int Nl){
    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);

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

    auto gr_DipMagZ = new TGraph();
    auto gr_DipMagZ_dd = new TGraph();

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

        gr_DipMagZ->SetPoint(i-Nl_0,time[i],DipMagZ[i]);
        gr_DipMagZ_dd->SetPoint(i-Nl_0,time[i],DipMagZ_dd[i]);
    }
    gr_CurS->SetTitle(";Time;Current at x=0");
    gr_DipX->SetTitle(";Time;p_{x}");
    gr_DipY->SetTitle(";Time;p_{y}");
    gr_DipX_dd->SetTitle(";Time;d^{2}p_{x}/dt^{2}");
    gr_DipY_dd->SetTitle(";Time;d^{2}p_{y}/dt^{2}");
    gr_QuadXX->SetTitle(";Time;Q_{xx}");
    gr_QuadXY->SetTitle(";Time;Q_{xy}");
    gr_QuadYY->SetTitle(";Time;Q_{yy}");
    gr_QuadXX_ddd->SetTitle(";Time;d^{3}Q_{xx}/dt^{3}");
    gr_QuadXY_ddd->SetTitle(";Time;d^{3}Q_{xy}/dt^{3}");
    gr_QuadYY_ddd->SetTitle(";Time;d^{3}Q_{yy}/dt^{3}");
    gr_DipMagZ->SetTitle(";Time;m_{z}");
    gr_DipMagZ_dd->SetTitle(";Time;d^{2}m_{z}/dt^{2}");

    c1->cd();   gr_CurS->Draw();c1->SaveAs("./CurS_t.pdf");

    c1->Clear(); gr_DipX->Draw();c1->SaveAs("./DipX_t.pdf");
    c1->Clear(); gr_DipY->Draw();c1->SaveAs("./DipY_t.pdf");
    c1->Clear(); gr_DipX_dd->Draw();c1->SaveAs("./DipX_dd_t.pdf");
    c1->Clear(); gr_DipY_dd->Draw();c1->SaveAs("./DipY_dd_t.pdf");

    c1->Clear(); gr_QuadXX->Draw();c1->SaveAs("./QuadXX_t.pdf");
    c1->Clear(); gr_QuadXY->Draw();c1->SaveAs("./QuadXY_t.pdf");
    c1->Clear(); gr_QuadYY->Draw();c1->SaveAs("./QuadYY_t.pdf");
    c1->Clear(); gr_QuadXX_ddd->Draw();c1->SaveAs("./QuadXX_ddd_t.pdf");
    c1->Clear(); gr_QuadXY_ddd->Draw();c1->SaveAs("./QuadXY_ddd_t.pdf");
    c1->Clear(); gr_QuadYY_ddd->Draw();c1->SaveAs("./QuadYY_ddd_t.pdf");

    c1->Clear(); gr_DipMagZ->Draw();c1->SaveAs("./DipMag_t.pdf");
    c1->Clear(); gr_DipMagZ_dd->Draw();c1->SaveAs("./DipMagZ_dd_t.pdf");

    delete gr_CurS; delete gr_DipX; delete gr_DipY; delete gr_DipX_dd; delete gr_DipY_dd;
    delete gr_QuadXX; delete gr_QuadXY; delete gr_QuadYY; delete gr_QuadXX_ddd; delete gr_QuadXY_ddd; delete gr_QuadYY_ddd;
    delete gr_DipMagZ; delete gr_DipMagZ_dd;
    delete c1;
}

void Radiation::PrintGraphInterpolation(int Nl_0, int Nl){
    auto c1 = new TCanvas("c1", "c1", 200, 10, 1350, 720);

    auto grI_CurS = new TGraph();

    auto grI_DipX_dd = new TGraph();
    auto grI_DipY_dd = new TGraph();

    auto grI_QuadXX_ddd = new TGraph();
    auto grI_QuadXY_ddd = new TGraph();
    auto grI_QuadYY_ddd = new TGraph();

    auto grI_DipMagZ_dd = new TGraph();

    int Np_interpolation = 10000;
    double t_interpolation = time[Nl_0];
    for(int i=0; i<Np_interpolation; i++){
        grI_DipX_dd->SetPoint(i,t_interpolation,SDipX_dd_t.Interpolate(t_interpolation));
        grI_DipY_dd->SetPoint(i,t_interpolation,SDipY_dd_t.Interpolate(t_interpolation));

        grI_QuadXX_ddd->SetPoint(i,t_interpolation,SQuadXX_ddd_t.Interpolate(t_interpolation));
        grI_QuadXY_ddd->SetPoint(i,t_interpolation,SQuadXY_ddd_t.Interpolate(t_interpolation));
        grI_QuadYY_ddd->SetPoint(i,t_interpolation,SQuadYY_ddd_t.Interpolate(t_interpolation));

        grI_DipMagZ_dd->SetPoint(i,t_interpolation,SDipMagZ_dd_t.Interpolate(t_interpolation));

        t_interpolation += (time[Nl]-time[Nl_0])/Np_interpolation;
    }
    grI_DipX_dd->SetTitle(";Time;d^{2}p_{x}/dt^{2}");
    grI_DipY_dd->SetTitle(";Time;d^{2}p_{y}/dt^{2}");
    grI_QuadXX_ddd->SetTitle(";Time;d^{3}Q_{xx}/dt^{3}");
    grI_QuadXY_ddd->SetTitle(";Time;d^{3}Q_{xy}/dt^{3}");
    grI_QuadYY_ddd->SetTitle(";Time;d^{3}Q_{yy}/dt^{3}");
    grI_DipMagZ_dd->SetTitle(";Time;d^{2}m_{z}/dt^{2}");

    c1->Clear(); grI_DipX_dd->Draw("AP");c1->SaveAs("./IDipX_dd_t.pdf");
    c1->Clear(); grI_DipY_dd->Draw("AP");c1->SaveAs("./IDipY_dd_t.pdf");

    c1->Clear(); grI_QuadXX_ddd->Draw("AP");c1->SaveAs("./IQuadXX_ddd_t.pdf");
    c1->Clear(); grI_QuadXY_ddd->Draw("AP");c1->SaveAs("./IQuadXY_ddd_t.pdf");
    c1->Clear(); grI_QuadYY_ddd->Draw("AP");c1->SaveAs("./IQuadYY_ddd_t.pdf");

    c1->Clear(); grI_DipMagZ_dd->Draw("AP");c1->SaveAs("./IDipMagZ_dd_t.pdf");

    delete grI_CurS; delete grI_DipX_dd; delete grI_DipY_dd; delete grI_QuadXX_ddd; delete grI_QuadXY_ddd; delete grI_QuadYY_ddd; delete grI_DipMagZ_dd;

    delete c1;
}


//Fields//
Vec Radiation::E_field(double t, double x, double y, double z){
    Vec E_Dip = Vec(3,0.); Vec E_Quad = Vec(3,0.); Vec E_DipMag = Vec(3,0.);

    double theta_i = theta_i_original;

    double r2 = x*x+y*y+z*z;
    double r = sqrt(r2);
    double t_r = t-r/(300*vF);
    if(t_r < t_min || t_r > t_max) cout << "Working outsided interpolated points [t_retarded = " << t_r << " ]\n";

    double* array = new double[3]{(x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*x/r2-SDipX_dd_t.Interpolate(t_r),
                      (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*y/r2-SDipY_dd_t.Interpolate(t_r),
                      (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*z/r2};
    E_Dip.SetEntries(3,array);
    delete array;

    double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
    double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
    double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

    double* array2 = new double[3]{(x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X,
                      (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y,
                      (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z};
    E_Quad.SetEntries(3,array2);
    delete array2;

    double* array3 = new double[3]{-y*SDipMagZ_dd_t.Interpolate(t_r), x*SDipMagZ_dd_t.Interpolate(t_r), 0};
    E_DipMag.SetEntries(3,array3);
    delete array3;

    return (1/(4*M_PI*r))*E_Dip*rotation180 + (1/(24*M_PI*300*vF*r))*E_Quad + (1/(4*M_PI*r2*300*vF))*E_DipMag*rotation180;
}

Vec Radiation::E_field_Image(double t, double x, double y, double z){
    Vec E_Dip = Vec(3,0.); Vec E_Quad = Vec(3,0.); Vec E_DipMag = Vec(3,0.); Vec E_field_v = Vec(3,0.); Vec E_field_T = Vec(3,0.);

    double theta_r = theta_r_image; 

    double r2 = x*x+y*y+z*z;
    double r = sqrt(r2);

    double theta_i = asin(sin(theta_r)/n_index);

    double r_vacuum;
    double r_dieletric = 2*D_Conductor/cos(theta_i);
    if(theta_r>1.5707) r_vacuum = d_xy;
    else r_vacuum = z/cos(theta_r);
    double t_r = t - (n_index*r_dieletric+r_vacuum)/(300*vF);
    if(t_r < t_min || t_r > t_max) cout << "Working outsided interpolated points [t_retarded = " << t_r << " ]\n";

    double* array = new double[3]{(x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*x/r2-SDipX_dd_t.Interpolate(t_r),
                      (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*y/r2-SDipY_dd_t.Interpolate(t_r),
                      (x*SDipX_dd_t.Interpolate(t_r)+y*SDipY_dd_t.Interpolate(t_r))*z/r2};
    E_Dip.SetEntries(3,array);
    delete array;

    double Q_vec_X = (SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r;
    double Q_vec_Y = (SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r;
    double Q_vec_Z = -(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r;

    double* array2 = new double[3]{(x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*x/r2-Q_vec_X,
                      (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*y/r2-Q_vec_Y,
                      (x*Q_vec_X+y*Q_vec_Y+z*Q_vec_Z)*z/r2-Q_vec_Z};
    E_Quad.SetEntries(3,array2);
    delete array2;

    double* array3 = new double[3]{-y*SDipMagZ_dd_t.Interpolate(t_r), x*SDipMagZ_dd_t.Interpolate(t_r), 0};
    E_DipMag.SetEntries(3,array3);
    delete array3;

    E_field_v = (1/(4*M_PI*r))*E_Dip*rotation180*mu_r + (1/(24*M_PI*300*vF*r))*E_Quad*mu_r*n_index + (1/(4*M_PI*r2*300*vF))*E_DipMag*rotation180*mu_r*n_index;

    double t_T, t_ll;
    if(theta_i<1e-3){
        t_T  = 2*n_index/(mu_r+n_index);
        t_ll = t_T;
    }else{
        t_T  = sqrt(cos(theta_r)/cos(theta_i))*2*n_index*cos(theta_i)/(mu_r*cos(theta_i)+n_index*cos(theta_r)); //perpendicular transmited coeficient
        t_ll = sqrt(cos(theta_r)/cos(theta_i))*2*n_index*cos(theta_i)/(mu_r*cos(theta_r)+n_index*cos(theta_i)); //parallel transmited coeficient
    }

    double* array_T = new double[3]{E_field_v[0]*sin(theta), E_field_v[1]*cos(theta), 0};
    E_field_T.SetEntries(3,array_T);
    delete array_T;

    return -(E_field_T*t_T + (E_field_v-E_field_T)*t_ll);
}

Vec Radiation::H_field(double t, double x, double y, double z){
    Vec r_vec = Vec(3,0.); Vec p_dd = Vec(3,0.); Vec Q_ddd = Vec(3,0.); Vec H_DipMag = Vec(3,0.);

    double r2 = x*x+y*y+z*z;
    double r = sqrt(r2);
    double t_r = t-r/(300*vF);

    double* array_r = new double[3]{x,y,z};
    r_vec.SetEntries(3,array_r);
    delete array_r;

    double* array_p = new double[3]{SDipX_dd_t.Interpolate(t_r),SDipY_dd_t.Interpolate(t_r),0};
    p_dd.SetEntries(3,array_p);
    delete array_p;

    double* array_Q = new double[3]{(SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r
                        ,(SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r
                        ,-(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r};
    Q_ddd.SetEntries(3,array_Q);
    delete array_Q;

    double* array3 = new double[3]{z*x*SDipMagZ_dd_t.Interpolate(t_r)/r2
                        , z*y*SDipMagZ_dd_t.Interpolate(t_r)/r2
                        , SDipMagZ_dd_t.Interpolate(t_r)*(z*z/r2-1)};
    H_DipMag.SetEntries(3,array3);
    delete array3;

    return (-1/(4*M_PI*r2))*r_vec.ex(p_dd)*rotation180 + (-1/(24*M_PI*300*vF*r2))*r_vec.ex(Q_ddd) + (1/(4*M_PI*r*300*vF))*H_DipMag*rotation180;
}

Vec Radiation::H_field_Image(double t, double x, double y, double z){
    Vec r_vec = Vec(3,0.); Vec p_dd = Vec(3,0.); Vec Q_ddd = Vec(3,0.); Vec H_DipMag = Vec(3,0.);
    Vec H_field_v = Vec(3,0.); Vec H_field_T = Vec(3,0.);

    double theta_r = theta_r_image; 

    double r2 = x*x+y*y+z*z;
    double r = sqrt(r2);

    double theta_i = asin(sin(theta_r)/n_index);

    double r_vacuum;
    double r_dieletric = 2*D_Conductor/cos(theta_i);
    if(theta_r>1.5707) r_vacuum = d_xy;
    else r_vacuum = z/cos(theta_r);
    double t_r = t - (n_index*r_dieletric+r_vacuum)/(300*vF);

    double* array_r = new double[3]{x,y,z};
    r_vec.SetEntries(3,array_r);
    delete array_r;

    double* array_p = new double[3]{SDipX_dd_t.Interpolate(t_r),SDipY_dd_t.Interpolate(t_r),0};
    p_dd.SetEntries(3,array_p);
    delete array_p;

    double* array_Q = new double[3]{(SQuadXX_ddd_t.Interpolate(t_r)*x + SQuadXY_ddd_t.Interpolate(t_r)*y)/r
                        ,(SQuadXY_ddd_t.Interpolate(t_r)*x + SQuadYY_ddd_t.Interpolate(t_r)*y)/r
                        ,-(SQuadXX_ddd_t.Interpolate(t_r) + SQuadYY_ddd_t.Interpolate(t_r))*z/r};
    Q_ddd.SetEntries(3,array_Q);
    delete array_Q;

    double* array3 = new double[3]{z*x*SDipMagZ_dd_t.Interpolate(t_r)/r2
                        , z*y*SDipMagZ_dd_t.Interpolate(t_r)/r2
                        , SDipMagZ_dd_t.Interpolate(t_r)*(z*z/r2-1)};
    H_DipMag.SetEntries(3,array3);
    delete array3;

    H_field_v = (-1/(4*M_PI*r2))*r_vec.ex(p_dd)*rotation180*n_index + (-1/(24*M_PI*300*vF*r2))*r_vec.ex(Q_ddd)*n_index*n_index + (1/(4*M_PI*r*300*vF))*H_DipMag*rotation180*n_index*n_index;

    double t_T, t_ll;
    if(theta_i<1e-3){
        t_T  = 2*n_index/(mu_r+n_index);
        t_ll = t_T;
    }else{
        t_T  = sqrt(cos(theta_r)/cos(theta_i))*2*n_index*cos(theta_i)/(mu_r*cos(theta_i)+n_index*cos(theta_r)); //perpendicular transmited coeficient
        t_ll = sqrt(cos(theta_r)/cos(theta_i))*2*n_index*cos(theta_i)/(mu_r*cos(theta_r)+n_index*cos(theta_i)); //parallel transmited coeficient
    }

    double* array_T = new double[3]{H_field_v[0]*sin(theta), H_field_v[1]*cos(theta), 0};
    H_field_T.SetEntries(3,array_T);
    delete array_T;

    return -(H_field_T*t_T + (H_field_v-H_field_T)*t_ll);
}

pair<Vec,Vec> Radiation::Emitter(double t, double O_x, double O_y, double x, double y, double z){
    d_xy = sqrt((x-O_x)*(x-O_x)+(y-O_y)*(y-O_y));
    global_z = z;

    double d_star;

    if(n_index < 1+1e-5) n_index = 1+1e-5;

    if (z<0) return make_pair((*this).E_field(t,x-O_x,y-O_y,z),(*this).H_field(t,x-O_x,y-O_y,z));

    if(fabs(z) < 20){
        theta_r_image = M_PI/2;
        d_star = D_Conductor*tan(asin(sin(theta_r_image)/n_index));
    }else if(d_xy<10){
        theta_r_image = 0;
        d_star = D_Conductor*tan(asin(sin(theta_r_image)/n_index));
    }else{
        theta_r_image = NewtonsMethod();
        d_star = D_Conductor*tan(asin(sin(theta_r_image)/n_index));
    }

    if(fabs(x) < 1e-3)
        theta = M_PI/2;
    else 
        theta = fabs(atan((y - O_y)/(x - O_x)));

    if((InConductor(d_star*cos(theta)-O_x, d_star*sin(theta)-O_y) == 0) || z < 0){
        return make_pair((*this).E_field(t,x-O_x,y-O_y,z),(*this).H_field(t,x-O_x,y-O_y,z));
    }else{
        return make_pair((*this).E_field(t,x-O_x,y-O_y,z)+(*this).E_field_Image(t,x-O_x,y-O_y,z),(*this).H_field(t,x-O_x,y-O_y,z)+(*this).H_field_Image(t,x-O_x,y-O_y,z));
    }
}

pair<Vec,Vec> Radiation::EmittingGrid(double t, vector<double> O_x, vector<double> O_y, double x, double y, double z){
    Vec Vec_temp = Vec(3,0.);
    int number_of_emitters = O_x.size();
    if(number_of_emitters != O_y.size()){
        cout << "Position for emitters is not valid [insert 2 vectors where the first is x and the second y]" << endl;
        Vec temp;
        return make_pair(temp, temp);
    }
    pair<Vec,Vec> fields;
    pair<Vec,Vec> fields_temp;

    fields.first = Vec_temp;
    fields.second = Vec_temp;
    for(int i=0; i<number_of_emitters; i++){
        fields_temp = (*this).Emitter(t,O_x[i],O_y[i],x,y,z);
        fields.first += fields_temp.first;
        fields.second += fields_temp.second;
    }

    return make_pair(fields.first,fields.second);
}

pair<Vec,Vec> Radiation::EmittingGrid180Rotation(double t, vector<double> O_x, vector<double> O_y, vector<bool> O_rot, double x, double y, double z){
    Vec Vec_temp = Vec(3,0.);
    int number_of_emitters = O_x.size();
    if(number_of_emitters != O_y.size()){
        cout << "Position for emitters is not valid [insert 2 vectors where the first is x and the second y]" << endl;
        Vec temp;
        return make_pair(temp, temp);
    }
    pair<Vec,Vec> fields;
    pair<Vec,Vec> fields_temp;

    fields.first = Vec_temp;
    fields.second = Vec_temp;
    for(int i=0; i<number_of_emitters; i++){
        if(O_rot[i]==1) rotation180 = -1;
        else rotation180 = 1;
        fields_temp = (*this).Emitter(t,O_x[i],O_y[i],x,y,z);
        fields.first += fields_temp.first;
        fields.second += fields_temp.second;
    }
    rotation180 = 1;

    return make_pair(fields.first,fields.second);
}


//RadiationPattern//
void Radiation::RadiationPatternTimeAverage180(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);
    ofstream file;

    file.open ("S_180" + configuration + ".txt");
    file << "sphere of radius" << distance_measure << "\n";
    file << "x    y    z    |S|*distance^2\n";

    N = N*2;

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);
        if(cos(mu)<1e-5) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid(phase, O_x_v,O_y_v,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";

            Ncount++;
        }
    }
   file.close();
}

double Radiation::RadiationPatternTimeAverage360(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("S_360" + configuration + ".txt");
    file << "Last line as integral over surface" << " -> sphere of radius" << distance_measure << "\n";
    file << "x    y    z    |S|*distance^2\n";    

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            S_integral += integral/N_t_points;

            file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";
            
            Ncount++;
        }
    }
    file << S_integral*4*M_PI/Ncount;
    file.close();

    return S_integral*4*M_PI/Ncount;
}

void Radiation::RadiationPatternTimeAverage180Graph3D(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t, bool print){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);
    ofstream file;

    if(print==1){
        file.open ("S_180" + configuration + ".txt");
        file << "sphere of radius" << distance_measure << "\n";
        file << "x    y    z    |S|*distance^2\n";
    }

    N = N*2;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
    auto gr3d_sphere = new TGraph2D(N);

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);
        if(cos(mu)<1e-5) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid(phase, O_x_v,O_y_v,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            if(print==1) file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";
            gr3d_sphere->SetPoint(Ncount, x*integral/N_t_points, y*integral/N_t_points, z*integral/N_t_points);

            Ncount++;
        }
    }

    c1->Clear();
    gr3d_sphere->Draw("PCOL Fi");
    c1->SaveAs(("S_180" + configuration + ".pdf").c_str());

    if(print==1) file.close();
    delete gr3d_sphere, c1;
}

double Radiation::RadiationPatternTimeAverage360Graph3D(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t, bool print){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    if(print==1){
        file.open ("S_360" + configuration + ".txt");
        file << "Last line as integral over surface" << " -> sphere of radius" << distance_measure << "\n";
        file << "x    y    z    |S|*distance^2\n";
    }
    

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
    auto gr3d_4pi = new TGraph2D(N);

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            S_integral += integral/N_t_points;

            if(print==1) file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";
            gr3d_4pi->SetPoint(Ncount, x*integral/N_t_points, y*integral/N_t_points, z*integral/N_t_points);
            
            Ncount++;
        }
    }
    c1->Clear();
    gr3d_4pi->Draw("PCOL Fi");
    c1->SaveAs(("S_360" + configuration + ".pdf").c_str());

    if(print==1){
        file << S_integral*4*M_PI/Ncount;
        file.close();
    }

    delete gr3d_4pi, c1;

    return S_integral*4*M_PI/Ncount;
}

double Radiation::RadiationTimeAverageIntegralSmallAngle(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t, double small_angle){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    N = (int)(2.*N/(1-cos(small_angle)));

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        if(mu > small_angle) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }
            S_integral += integral/N_t_points;

            Ncount++;
        }
    }

    return 2*M_PI*(1-cos(small_angle))*S_integral/Ncount;
}

void Radiation::RadiationTimeAverageE_H_Planes180(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    double x1,y1,x2,y2,z,phi,S_mod1, S_mod2;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("E_H_Planes" + configuration + ".txt");
    file << "phi    |S|(E plane)    |S|(H plane)" << "\n";  

    Ncount = 0;

    phi = -M_PI/2;
    for(int i=0; i<N; i++){
        phi += 2*M_PI/N;

        x1 = sin(phi)*cos(E_plane_angle);
        y1 = sin(phi)*sin(E_plane_angle);
        x2 = sin(phi)*cos(E_plane_angle+M_PI/2);
        y2 = sin(phi)*sin(E_plane_angle+M_PI/2);
        z  = cos(phi);

        if(z<1e-4){
            continue; 
        } 

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x1,distance_measure*y1,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod1 = integral/N_t_points;
        Ncount++;

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x2,distance_measure*y2,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod2 = integral/N_t_points;

        file << phi << "   " << S_mod1 << "   " << S_mod2 << "\n";

        Ncount++;
    }
    file.close();
}

void Radiation::RadiationTimeAverageE_H_Planes360(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    double x1,y1,x2,y2,z,phi,S_mod1, S_mod2;
    double integral, phase, S_integral;

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("E_H_Planes" + configuration + ".txt");
    file << "phi    |S|(E plane)    |S|(H plane)" << "\n";

    Ncount = 0;

    phi = -M_PI;
    for(int i=0; i<N; i++){
        phi += 2*M_PI/N;

        x1 = sin(phi)*cos(E_plane_angle);
        y1 = sin(phi)*sin(E_plane_angle);
        x2 = sin(M_PI/2+M_PI/1000)*cos(phi);
        y2 = sin(M_PI/2+M_PI/1000)*sin(phi);
        z  = cos(phi);

        if(fabs(z)<1e-4){
            continue; 
        } 

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x1,distance_measure*y1,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod1 = integral/N_t_points;
        Ncount++;

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid(phase,O_x_v,O_y_v,distance_measure*x2,distance_measure*y2,distance_measure*cos(M_PI/2+M_PI/1000));

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod2 = integral/N_t_points;

        file << phi << "   " << S_mod1 << "   " << S_mod2 << "\n";

        Ncount++;
    }
    file.close();
}


void Radiation::RadiationPatternTimeAverage180_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);
    ofstream file;

    file.open ("S_180" + configuration + ".txt");
    file << "sphere of radius" << distance_measure << "\n";
    file << "x    y    z    |S|*distance^2\n";

    N = N*2;

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);
        if(cos(mu)<1e-5) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid180Rotation(phase, O_x_v,O_y_v,O_rot,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";

            Ncount++;
        }
    }
   file.close();
}

double Radiation::RadiationPatternTimeAverage360_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("S_360" + configuration + ".txt");
    file << "Last line as integral over surface" << " -> sphere of radius" << distance_measure << "\n";
    file << "x    y    z    |S|*distance^2\n";    

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            S_integral += integral/N_t_points;

            file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";

            Ncount++;
        }
    }
    file << S_integral*4*M_PI/Ncount;
    file.close();

    return S_integral*4*M_PI/Ncount;
}

void Radiation::RadiationPatternTimeAverage180Graph3D_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t, bool print){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);
    ofstream file;

    if(print==1){
        file.open ("S_180" + configuration + ".txt");
        file << "sphere of radius" << distance_measure << "\n";
        file << "x    y    z    |S|*distance^2\n";
    }

    N = N*2;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
    auto gr3d_sphere = new TGraph2D(N);

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);
        if(cos(mu)<1e-5) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid180Rotation(phase, O_x_v,O_y_v,O_rot,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            if(print==1) file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";
            gr3d_sphere->SetPoint(Ncount, x*integral/N_t_points, y*integral/N_t_points, z*integral/N_t_points);

            Ncount++;
        }
    }

    c1->Clear();
    gr3d_sphere->Draw("PCOL Fi");
    c1->SaveAs(("S_180" + configuration + ".pdf").c_str());

    if(print==1) file.close();
    delete gr3d_sphere, c1;
}

double Radiation::RadiationPatternTimeAverage360Graph3D_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t, bool print){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    if(print==1){
        file.open ("S_360" + configuration + ".txt");
        file << "Last line as integral over surface" << " -> sphere of radius" << distance_measure << "\n";
        file << "x    y    z    |S|*distance^2\n";
    }
    

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    auto c1 = new TCanvas("c1", "c1", 200, 10, 1280, 720);
    auto gr3d_4pi = new TGraph2D(N);

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }

            S_integral += integral/N_t_points;

            if(print==1) file << x*distance_measure << "   " << y*distance_measure << "   " << z*distance_measure << "   " << integral/N_t_points << "\n";
            gr3d_4pi->SetPoint(Ncount, x*integral/N_t_points, y*integral/N_t_points, z*integral/N_t_points);

            Ncount++;
        }
    }
    c1->Clear();
    gr3d_4pi->Draw("PCOL Fi");
    c1->SaveAs(("S_360" + configuration + ".pdf").c_str());

    if(print==1){
        file << S_integral*4*M_PI/Ncount;
        file.close();
    }

    delete gr3d_4pi, c1;

    return S_integral*4*M_PI/Ncount;
}

double Radiation::RadiationTimeAverageIntegralSmallAngle_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t, double small_angle){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    int M_phi, M_mu;
    double a ,d_sphere, d_mu, d_phi, mu, phi, x, y, z;
    double integral, phase, S_integral;

    N = (int)(2.*N/(1-cos(small_angle)));

    a = 4*M_PI/N;
    d_sphere = sqrt(a);
    M_mu = round(M_PI/d_sphere);
    d_mu = M_PI/M_mu;
    d_phi = a/d_mu;

    Ncount = 0;
    S_integral = 0;
    for(int m=0; m<M_mu-1; m++){
        mu = M_PI*(m + 0.5)/M_mu;
        M_phi = round(2*M_PI*sin(mu)/d_phi);

        if(mu > small_angle) continue;

        for(int n=0; n<M_phi-1; n++){
            phi = 2*M_PI*n/M_phi;
            x = sin(mu)*cos(phi);
            y = sin(mu)*sin(phi);
            z = cos(mu);

            integral = 0;
            for(int j=0; j<N_t_points; j++){
                phase = t0+(double)j*delta_t/N_t_points;

                fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x,distance_measure*y,distance_measure*z);

                if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
                else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            }
            S_integral += integral/N_t_points;

            Ncount++;
        }
    }

    return 2*M_PI*(1-cos(small_angle))*S_integral/Ncount;
}

void Radiation::RadiationTimeAverageE_H_Planes180_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    double x1,y1,x2,y2,z,phi,S_mod1, S_mod2;
    double integral, phase, S_integral;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("E_H_Planes" + configuration + ".txt");
    file << "phi    |S|(E plane)    |S|(H plane)" << "\n";  

    Ncount = 0;

    phi = -M_PI/2;
    for(int i=0; i<N; i++){
        phi += 2*M_PI/N;

        x1 = sin(phi)*cos(E_plane_angle);
        y1 = sin(phi)*sin(E_plane_angle);
        x2 = sin(phi)*cos(E_plane_angle+M_PI/2);
        y2 = sin(phi)*sin(E_plane_angle+M_PI/2);
        z  = cos(phi);

        if(z<1e-4){
            continue; 
        } 

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x1,distance_measure*y1,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod1 = integral/N_t_points;
        Ncount++;

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x2,distance_measure*y2,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod2 = integral/N_t_points;

        file << phi << "   " << S_mod1 << "   " << S_mod2 << "\n";

        Ncount++;
    }
    file.close();
}

void Radiation::RadiationTimeAverageE_H_Planes360_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t){
    pair<Vec,Vec> fieldsEH;
    int Ncount;
    double x1,y1,x2,y2,z,phi,S_mod1, S_mod2;
    double integral, phase, S_integral;

    if(Radius_Conductor>1e-3 || n_index> 1+1e-3) cout << "For a 360 degree calculation, refractive index needs to be 1 and no conductor (proceeding that way)\n";
    n_index = 1;
    Radius_Conductor = 0;

    ofstream file;
    string configuration = "_n=" + to_string_with_precision(n_index,1) + "_D_Cond=" + to_string_with_precision(D_Conductor,1) + "_R_Cond=" + to_string_with_precision(Radius_Conductor,1);

    file.open ("E_H_Planes" + configuration + ".txt");
    file << "phi    |S|(E plane)    |S|(H plane)" << "\n";  

    Ncount = 0;

    phi = -M_PI;
    for(int i=0; i<N; i++){
        phi += 2*M_PI/N;

        x1 = sin(phi)*cos(E_plane_angle);
        y1 = sin(phi)*sin(E_plane_angle);
        x2 = sin(phi)*cos(E_plane_angle+M_PI/2);
        y2 = sin(phi)*sin(E_plane_angle+M_PI/2);
        z  = cos(phi);

        if(fabs(z)<1e-4){
            continue; 
        } 

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x1,distance_measure*y1,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod1 = integral/N_t_points;
        Ncount++;

        integral = 0;
        for(int j=0; j<N_t_points; j++){
            phase = t0+(double)j*delta_t/N_t_points;

            fieldsEH = EmittingGrid180Rotation(phase,O_x_v,O_y_v,O_rot,distance_measure*x2,distance_measure*y2,distance_measure*z);

            if(j == 0 || j == N_t_points-1) integral += 0.5*distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
            else integral += distance_measure*distance_measure*Poyting(fieldsEH.first, fieldsEH.second).mod();
        }
        S_mod2 = integral/N_t_points;

        file << phi << "   " << S_mod1 << "   " << S_mod2 << "\n";

        Ncount++;
    }
    file.close();
}

//Sets//
void Radiation::SetVf(double x){
    if (vF<0) cout << "Fermi Velocity must be bigger than 0 \n";
    vF = x;
}
void Radiation::Set_n_index(double x){
    if (x < (1-1e-3)) cout << "Refractive index must be bigger than 1\n";
    n_index = x;
    mu_r = n_index*n_index/epsilon_r;
    if(mu_r < 1) cout << "relative permeability can't be smaller than 1 \n";
}
void Radiation::Set_epsilon_r(double x){
    if (x < (1-1e-3)) cout << "Refractive index must be bigger than 1\n";
    epsilon_r = x;
    mu_r = n_index*n_index/epsilon_r;
    if(mu_r < 1) cout << "relative permeability can't be smaller than 1 \n";
}
void Radiation::SetD_Conductor(double x){
    if (x<0) cout << "Distance to conductor must be bigger than 0 \n";
    D_Conductor = x;
}
void Radiation::SetR_Conductor(double x){
    if (x<0) cout << "Radius must not be smaller than 0 \n";
    Radius_Conductor = x;
}

void Radiation::SetDistance_Measure(double x){
    if (x<0) cout << "Measurement distance must be bigger than 0 \n";
    distance_measure = x;
}


//Gets//
double Radiation::GetVf(){
    return vF;
}
double Radiation::Get_n_index(){
    return n_index;
}
double Radiation::Get_epsilon_r(){
    return epsilon_r;
}
double Radiation::Get_mu_r(){
    return  mu_r;
}
double Radiation::GetD_Conductor(){
    return D_Conductor;
}
double Radiation::GetR_Conductor(){
    return Radius_Conductor;
}
double Radiation::GetDistance_Measure(){
    return distance_measure;
}


//Others//
double Radiation::NewtonsMethod(double x1, double x2, double error_min, int iter_max){
    int i = 0;
    double fx1 = 0., fx2 = 0., x3 = 0.;
    while((fabs(x1-x2) >= error_min || fabs(fx2) > error_min) && i < iter_max){
        fx2 = (*this).Theta_refraction_equation(x2);
        fx1 = (*this).Theta_refraction_equation(x1);
        if(fx2-fx1 != 0) x3 = x2 - (fx2*(x2-x1)/(fx2-fx1));
        else{
            cout << "We have found a stationary point at x = " << x2 << ". Ending Newton's method here." << endl;
            cout << "returning pi/2" << endl;
            return M_PI/2;
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

double Radiation::Theta_refraction_equation(double theta_r){
    return global_z*tan(theta_r) + 2*D_Conductor*sin(theta_r)/cos(asin(sin(theta_r)/n_index)) - d_xy;
}

bool Radiation::InConductor(double x, double y){
    if(x*x + y*y < (Radius_Conductor*Radius_Conductor))
        return 1;
    else
        return 0;
}


//foreign//
Vec Poyting(Vec E, Vec H){
    return E.ex(H);
}

std::string to_string_with_precision(const double a_value, const int n){
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}