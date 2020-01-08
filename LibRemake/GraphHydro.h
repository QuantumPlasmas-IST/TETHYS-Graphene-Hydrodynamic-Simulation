void BoundaryCond(int type, int N, float * den, float * vel); //
void ExtremaFinding(float * vec_in, int N, float sound, float dt,float & sat, float  & tau, float & error, std::string extremafile);
float ImagFreq(float sound, float fermi, float col_freq);                  //
void InitialCondSine(int N, float dx,  float * den, float * vel); //
void InitialCondRand(int N, float dx,  float * den, float * vel); //
float PhaseVel(float sound, float fermi);
float RealFreq(float sound, float fermi, float col_freq, int mode);  //
void ShockFinding(float * in, int N, float t , float dx,  std::string shockfile);
float TimeStepCFL(float dx, float sound, float fermi);
