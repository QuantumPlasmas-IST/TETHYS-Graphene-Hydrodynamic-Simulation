/************************************************************************************************\
* 2022 Rui Martins, Pedro Cosme                                                                 *
*                                                                                               *
*                                                                                               *
\************************************************************************************************/

#include "../src/includes/RadiationLib.h"

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

	Radiation Rad(10.5,2.5,1.5,2,0.05,3000);

	int Nl = Rad.ReadFromElectroFile(input_file_name);

    int Nl_0 = 1700;

    Rad.Interpolation(Nl_0,Nl);

    //Rad.PrintGraphElectro(Nl_0, Nl);
    //Rad.PrintGraphInterpolation(Nl_0, Nl);

    vector<double> O_x1{0};
    vector<double> O_y1{0};

    vector<double> O_x2{0,1};
    vector<double> O_y2{0,0};

    vector<double> O_x5{0, 1, 0, -1, 0};
    vector<double> O_y5{0, 0, 1, 0, -1};

    vector<double> O_xN, O_yN;
    for(int i=-8; i<=8; i++){
        for(int j=-8; j<=8; j++){
            O_xN.push_back((double)i);
            O_yN.push_back((double)j);
        }
    }

    /*/
    Rad.Set_n_index(2.5);
    Rad.SetD_Conductor(30);
	Rad.SetR_Conductor(30);
    Rad.RadiationPatternTimeAverage180Graph(8.5, O_x1, O_y1, 50000, 100, 0.356, 1);
    //*/

    /*/
    Rad.Set_n_index(1.);
	Rad.SetD_Conductor(0);
	Rad.SetR_Conductor(0.);
    cout << Rad.RadiationPatternTimeAverage360Graph(8.5, O_x1, O_y1, 50000, 100, 0.356, 1) << endl;
    //*/
	
    Rad.Set_n_index(2.5);
	Rad.SetD_Conductor(1.5);
	Rad.SetR_Conductor(2);
	cout << Rad.RadiationTimeAverageIntegralSmallAngle(8.5, O_x5, O_y5, 100000, 100, 0.356, M_PI/18) << endl;
    
    cout << "[1A\033[2K\033[1;32mDONE!\033[0m\n";
    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
    cout << "═══════════════════════════════════════════════════════════════════════════" <<endl;
	return 0;
}