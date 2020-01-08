float DtElectricDipole(int N,float dx, float * cur); //
float TotalCurrent(int N,float dx, float * den,float * vel); //
float TotalElectricDipole(int N,float dx, float * den);      //
float KineticEnergy(int N,float dx, float * den, float * vel); //
float CartDistance(float x , float y , float z, float X , float Y , float Z );
float RetardedTime(float time, float x , float y , float z, float X , float Y , float Z );	 
void  JefimenkoEMField(int XDIM, int YDIM, float dx, float dy, float dt, float Xpos, float Ypos, float Zpos,  float ** rho, float ** rho_dot, float ** cur, float ** cur_dot, float Time , float  * E_out , float  * B_out, float  * S_out   );
