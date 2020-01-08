// HYDRO 
float TimeStepCFL(float dx, float sound, float fermi);
float PhaseVel(float sound, float fermi);
float RealFreq(float sound, float fermi, float col_freq, int mode);  //
float ImagFreq(float sound, float fermi, float col_freq);                  //
void BoundaryCond(int type, int N, float * den, float * vel); //
void InitialCondSine(int N, float dx,  float * den, float * vel); //
void InitialCondRand(int N, float dx,  float * den, float * vel); //
void ExtremaFinding(float * vec_in, int N, float sound, float dt,float & sat, float  & tau, float & error, std::string extremafile);
void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile);

// ELECTRO
float DtElectricDipole(int N,float dx, float * cur); //
float TotalCurrent(int N,float dx, float * den,float * vel); //
float TotalElectricDipole(int N,float dx, float * den);      //
float KineticEnergy(int N,float dx, float * den, float * vel); //
float CartDistance(float x , float y , float z, float X , float Y , float Z );
float RetardedTime(float time, float x , float y , float z, float X , float Y , float Z );	 
void  JefimenkoEMField(int XDIM, int YDIM, float dx, float dy, float dt, float Xpos, float Ypos, float Zpos,  float ** rho, float ** rho_dot, float ** cur, float ** cur_dot, float Time , float  * E_out , float  * B_out, float  * S_out   );

// GENERAL
float RootMeanSquare(int N, float dt, float * f);
float SignalAverage(int N, float dt, float * f);
float GaussKernel(int position , float t); //
float GaussKernelDerivative(int position , float t); //
void ConvolveGauss(int type, float M, float t, float * in, float * out, int size);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );
void Autocorrelation(float * out_gamma ,float * in , int crop, int size);
void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );

void BannerDisplay(void);
void WellcomeScreen(float vel_snd, float vel_fer,float col_freq, float dt, float dx, float Tmax);
void RecordLogFile(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float Tmax);
