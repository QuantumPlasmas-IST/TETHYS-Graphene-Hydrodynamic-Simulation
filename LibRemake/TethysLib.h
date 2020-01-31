void Autocorrelation(float * out_gamma ,float * in , int crop, int size);
void AverageFilter(float * vec_in, float * vec_out, int size , int width );
void BannerDisplay(void);
void ConvolveGauss(int type, float M, float t, float * in, float * out, int size);
float GaussKernel(int position , float t); //
float GaussKernelDerivative(int position , float t); //
void RecordLogFile(float vel_snd, float vel_fer, float col_freq, float dt, float dx, float Tmax);
float RootMeanSquare(int N, float dt, float * f);
float SignalAverage(int N, float dt, float * f);
void SpaceDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void TimeDerivative(int size_rows,int size_cols, float dt,float ** f_in , float ** df_out );
void WellcomeScreen(float vel_snd, float vel_fer,float col_freq, float dt, float dx, float Tmax);


void BoundaryCond(int type, int N, float * den, float * vel); //
void ExtremaFinding(float * vec_in, int N, float sound, float dt,float & sat, float  & tau, float & error, std::string extremafile);
float ImagFreq(float sound, float fermi, float col_freq);                  //
void InitialCondSine(int N, float dx,  float * den, float * vel); //
void InitialCondRand(int N, float dx,  float * den, float * vel); //
float PhaseVel(float sound, float fermi);
float RealFreq(float sound, float fermi, float col_freq, int mode);  //
void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile);
float TimeStepCFL(float dx, float sound, float fermi);



float DensityFlux(float den,float vel,float vel_snd,float vel_fer);

float VelocityFlux(float den,float vel,float vel_snd,float vel_fer);

float EnergyFlux(float den,float vel,float vel_snd,float vel_fer);

float DensitySource(float den,float vel,float vel_snd,float vel_fer);

float VelocitySource(float den,float vel,float vel_snd,float vel_fer,float col_freq);

float EnergySource(float den,float den_der,float vel,float vel_snd,float vel_fer);
