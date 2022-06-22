/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281                                                                     *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#ifndef RADIATIONLIB_H
#define RADIATIONLIB_H

#include <time.h>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "../Library_Radiation/Spline3Interpolator.h"

using namespace std;


class Radiation{
public:
    Radiation(double n_index_ref=2.5, double Dist_Conductor=1.5, double R_Conductor=2, double Dist_dieletric=0.05, double Distance_measure=3000);
    
    //Read a file comming from EletronicAnalysis.cpp
    int ReadFromElectroFile(const string input_file_name);
   
    //Cubic spline interpolation of the points stored and only on the used quatites
    void Interpolation(int Nl_0, int Nl); 

    Vec E_field(double t, double x, double y, double z);
    Vec E_field_Image(double t, double x, double y, double z);
    Vec H_field(double t, double x, double y, double z);
    Vec H_field_Image(double t, double x, double y, double z);
    pair<Vec,Vec> Emitter(double t, double O_x, double O_y, double x, double y, double z);
    pair<Vec,Vec> EmittingGrid(double t, vector<double> O_x, vector<double> O_y, double x, double y, double z);
    pair<Vec,Vec> EmittingGrid180Rotation(double t, vector<double> O_x, vector<double> O_y, vector<bool> O_rot, double x, double y, double z);

    //vectors O_x_v and O_y_v contain the position of each emmiter in the z=0 plane

    //creates a file with the positions on a semi-sphere of radius distance_measure and |Poyting|*radius^2
    void RadiationPatternTimeAverage180(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t);
    
    //returns total time averaged flux of Poyinting vector and creates a file with the positions on a sphere of radius distance_measure and |Poyting|*radius^2
    double RadiationPatternTimeAverage360(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t);

    //returns time averaged flux of Poyinting vector across the top angles
    double RadiationTimeAverageIntegralSmallAngle(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t, double small_angle);
    
    //creates a file with |Poyting|*radius^2 on the intersection of the semi-sphere of radius distance_measure with the planes "E" and "H"
    void RadiationTimeAverageE_H_Planes180(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t);
    
    //creates a file with |Poyting|*radius^2 on the intersection of the sphere of radius distance_measure with the planes "E" and "H"
    void RadiationTimeAverageE_H_Planes360(double t0, vector<double> O_x_v, vector<double> O_y_v, int N, int N_t_points, double delta_t);

    //Same as the above functions but allows for 180 rotations on the emitters
    //If in O_rot it's 1 emmiter will be considered rotated
    void RadiationPatternTimeAverage180_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t);
    double RadiationPatternTimeAverage360_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t);

    double RadiationTimeAverageIntegralSmallAngle_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t, double small_angle);
    void RadiationTimeAverageE_H_Planes180_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t);
    void RadiationTimeAverageE_H_Planes360_180Rotation(double t0, vector<double> O_x_v, vector<double> O_y_v, vector<bool> O_rot, int N, int N_t_points, double delta_t);


    void SetVf(double x);
    void Set_n_index(double x);
    void SetD_Conductor(double x);
    void SetR_Conductor(double x);
    void SetD_Dieletric(double x);
    void SetDistance_Measure(double x);

    double GetVf(double x);
    double Get_n_index(double x);
    double GetD_Conductor(double x);
    double GetR_Conductor(double x);
    double GetD_Dieletric(double x);
    double GetDistance_Measure(double x);

    double NewtonsMethod(double x1=1.56, double x2=1.57, double error_min=1e-5, int iter_max=40);
    double Theta_refraction_equation(double theta_r);

    bool InConductor(double x, double y);

private:
    string input_file_name;
    
    vector<double> time, CurS;
    vector<double> DipX, DipX_dd, DipY, DipY_dd;
    vector<double> DipMagZ, DipMagZ_dd; 
    vector<double> QuadXX, QuadXX_ddd, QuadXY, QuadXY_ddd, QuadYY, QuadYY_ddd;

    int rotation180 = 1;

    double vF, n_index;
    double D_Conductor, Radius_Conductor, D_dieletric, distance_measure;
    double d_xy, global_z;
    double theta, theta_i_original, theta_r_image;
    double E_plane_angle;

    Spline3Interpolator SDipX_dd_t = Spline3Interpolator();
    Spline3Interpolator SDipY_dd_t = Spline3Interpolator();
    Spline3Interpolator SDipMagZ_dd_t = Spline3Interpolator();
    Spline3Interpolator SQuadXX_ddd_t = Spline3Interpolator();
    Spline3Interpolator SQuadXY_ddd_t = Spline3Interpolator();
    Spline3Interpolator SQuadYY_ddd_t = Spline3Interpolator();
};

Vec Poyting(Vec E, Vec H);

string to_string_with_precision(const double a_value, const int n);

#endif
