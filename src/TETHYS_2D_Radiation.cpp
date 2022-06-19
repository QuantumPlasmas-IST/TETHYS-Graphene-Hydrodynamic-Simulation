/************************************************************************************************\
* 2022 Rui Martins, Pedro Cosme                                                                 *
*                                                                                               *
*                                                                                               *
\************************************************************************************************/

#include "includes/RadiationLib.h"

using namespace std;

int main(int argc, char **argv){
    string input_file_name;
    if(argc==2){
        input_file_name = argv[1];
    }else{
        cout <<"Input whole name of 'electro_2D_[...].dat' file to read:\t";
        cin >> input_file_name;
    }
    clock_t tStart = clock();

    //v_F is coming from the name of the file, if that is change use Rad.SetVf(double) here

    Radiation Rad(2.5,1.5,2,0.05,3000);

    //Returns last line of file
    int Nl = Rad.ReadFromElectroFile(input_file_name);

    int Nl_0 = 0;

    //We advise to interpolate only the interval that will be used to save computer power
    Rad.Interpolation(Nl_0,Nl);


    //array with only 1 emitter
    vector<double> O_x1{0};
    vector<double> O_y1{0};

    //array with 5x5 square grid of emitters and even rows rotated (using _180rotation functions)
    vector<double> O_xN, O_yN;
    vector<bool> O_rotN;
    for(int i=-2; i<=2; i++){
        for(int j=-2; j<=2; j++){
            O_xN.push_back((double)i);
            O_yN.push_back((double)j);
            if(j%2 == 1) O_rotN.push_back(1);
            else O_rotN.push_back(0);
        }
    }

    //See headers for more comments on functions
    //Using these functions the program will calculate EH fields (#time_points)*(#space_points)*(#emitters) times
    //Each minute should make 5e6 field calculations

    //Introducing 9.1 in time, in the retarded time will be ~ 9.1-distance_measure/(v_F*300)
    //It's the retarded time + time interval that needs to be inside of the interpolated interval


    //
    Rad.Set_n_index(1.0);
    Rad.SetD_Conductor(0.75);
    Rad.SetR_Conductor(0.0);
    cout << Rad.RadiationTimeAverageIntegralSmallAngle(9.1, O_x1, O_y1, 800, 100, 0.36, M_PI/18) << endl;
    cout << Rad.RadiationPatternTimeAverage360(9.1, O_x1, O_y1, 50000, 100, 0.36) << endl;
    Rad.RadiationTimeAverageE_H_Planes360(9.1, O_x1, O_y1, 200, 100, 0.36);
    //*/

    //
    Rad.Set_n_index(2.5);
    Rad.SetD_Conductor(0.5);
    Rad.SetR_Conductor(5);
    cout << Rad.RadiationTimeAverageIntegralSmallAngle_180Rotation(9.1, O_xN, O_yN, O_rotN, 800, 100, 0.36, M_PI/18) << "\n";
    Rad.RadiationPatternTimeAverage180_180Rotation(9.1, O_xN, O_yN, O_rotN, 3000, 100, 0.36);
    Rad.RadiationTimeAverageE_H_Planes180_180Rotation(9.1, O_xN, O_yN, O_rotN, 200, 100, 0.36);
    //*/
    

    cout << "[1A\033[2K\033[1;32mDONE!\033[0m\n";
    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
    return 0;
}