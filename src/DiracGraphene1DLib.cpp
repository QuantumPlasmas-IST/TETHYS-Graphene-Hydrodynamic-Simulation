

#include "includes/DiracGraphene1DLib.h"

DiracGraphene1D::DiracGraphene1D(SetUpParametersCNP &input_parameters) : Fluid1D(input_parameters) {

    hsize_t dimsf[1];
    dimsf[0] = static_cast<hsize_t>(Nx);
    DataspaceHDen = new DataSpace(RANK, dimsf );
    DataspaceHVel = new DataSpace(RANK, dimsf );
    GrpHDen = nullptr;
    GrpHVel = nullptr;

    vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
    col_freq = input_parameters.CollisionFrequency ; // collision_frequency
    cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency
    therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity

    vel_therm = input_parameters.ThermalVelocity;
    A = input_parameters.Diffusive_sourceterm;
    B = input_parameters.Creation_sourceterm;

    //param = {vel_snd,vel_fer,0.0f,kin_vis,0.0f,therm_diff,col_freq,0.0f};

    char buffer [100];
    sprintf (buffer, "S=%.2fvF=%.2fvT=%.2fA=%.3fB=%.3fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, vel_therm, A, B, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
    file_infix = buffer;

    // main grid variables Nx
    HDen = new float[Nx];
    HVel = new float[Nx];

    HoleUmain = new StateVec1D[Nx]();
    HoleUmid = new StateVec1D[Nx-1]();
}

/*
void DiracGraphene1D::SetSimulationTime() {
    this->SetTmax(5.0f+0.02f*vel_therm+20.0f/vel_therm);
}
*/

void DiracGraphene1D::CflCondition(){ // Eventual redefinition

    dx = lengX / ( float ) ( Nx - 1 );
    float lambda;
    if(vel_snd<0.36f*vel_fer){
        lambda=1.2f*vel_fer;
    }else{
        lambda=1.97f*vel_snd + 0.5f*vel_fer;
    }
    dt = dx/lambda;

    /*
    dx = lengX / ( float ) ( Nx - 1 );
    float lambda;
    if(vel_therm<vel_snd){
        lambda=0.75f+sqrt(1.6f+3.0f*vel_snd*vel_snd);
    }else{
        lambda=0.75f+sqrt(1.6f+3.0f*vel_therm*vel_therm);
    }
    float dt1 = dx/abs(lambda);     // retirei fator 0.5f, era suposto estar la?
    dt = dt1;


    //  CFL condition for (1,9) Weighted explicit method

    if(kin_vis>0.0){
        float dt2 = 0.4f*dx*dx/max(therm_diff,kin_vis);  // faz sentido ser o max entre os dois parametros?
        dt = min(dt1, dt2);                              // fonte: TETHYS formula (8)
    }
*/

}


void DiracGraphene1D::SetSound(){
    for(int i = 0; i<Nx  ;i++){
        vel_snd_arr[i]= vel_snd;
        Umain[i].S()=vel_snd;
        HoleUmain[i].S()=vel_snd;
        //Uaux[i].S()=vel_snd;
    }
    for(int i = 0; i<Nx-1  ;i++){
        Umid[i].S()=vel_snd;
        HoleUmid[i].S()=vel_snd;
    }
}


void DiracGraphene1D::InitialCondRand(){
    random_device rd;
    float maxrand;
    maxrand = (float) random_device::max();
    for (int c = 0; c < Nx; c++ ){
        float noise;
        noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        //Den[c] = 1.0f + 0.05f * (noise - 0.5f);
        Umain[c].n()=1.0f + 0.05f * (noise - 0.5f);
        noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        //HDen[c] = 1.0f + 0.05f * (noise - 0.5f);
        HoleUmain[c].n()=1.0f + 0.05f * (noise - 0.5f);
    }
    this -> SetSound();
}


void DiracGraphene1D::RichtmyerStep1() {
    for( int i = 0; i <= Nx -2; i++){

        StateVec1D eUavg{};
        eUavg = 0.5f*(Umain[i+1] + Umain[i]);

        StateVec1D hUavg{};
        hUavg = 0.5f*(HoleUmain[i+1] + HoleUmain[i]);

        //____________________ Electrons _________________________

        Umid[i].n() = eUavg.n() - 0.5f*(dt/dx)*(EleDensityFlux(Umain[i+1]) - EleDensityFlux(Umain[i]))
                      + (0.5f*dt) * EleDensitySource(eUavg, hUavg) ;

        Umid[i].v() = eUavg.v() - 0.5f*(dt/dx)*(EleVelocityFlux(Umain[i+1]) - EleVelocityFlux(Umain[i]))
                      + (0.5f*dt) * EleVelocitySource(eUavg) ;

        //____________________ Holes _________________________

        HoleUmid[i].n() = hUavg.n() - 0.5f*(dt/dx)*(HolDensityFlux(HoleUmain[i+1]) - HolDensityFlux(HoleUmain[i]))
                      + (0.5f*dt) * HolDensitySource(eUavg, hUavg) ;

        HoleUmid[i].v() = hUavg.v() - 0.5f*(dt/dx)*(HolVelocityFlux(HoleUmain[i+1]) - HolVelocityFlux(HoleUmain[i]))
                      + (0.5f*dt) * HolVelocitySource(hUavg) ;

    }
}

void DiracGraphene1D::RichtmyerStep2() {
    for( int i = 1; i <= Nx -2; i++ ){

        StateVec1D eUold(Umain[i]);
        StateVec1D hUold(HoleUmain[i]);

        //____________________ Electrons _________________________

        Umain[i].n() = eUold.n() - (dt/dx) * (EleDensityFlux(Umid[i]) - EleDensityFlux(Umid[i-1]))
                       + dt * EleDensitySource(eUold,hUold);

        Umain[i].v() = eUold.v() - (dt/dx)*(EleVelocityFlux(Umid[i]) - EleVelocityFlux(Umid[i-1]))
                           + dt * EleVelocitySource(eUold);

        //____________________ Holes _________________________

        HoleUmain[i].n() = hUold.n() - (dt/dx) * (HolDensityFlux(HoleUmid[i]) - HolDensityFlux(HoleUmid[i-1]))
                           + dt * HolDensitySource(eUold, hUold);

        HoleUmain[i].v() = hUold.v() - (dt/dx)*(HolVelocityFlux(HoleUmid[i]) - HolVelocityFlux(HoleUmid[i-1]))
                                       + dt * HolVelocitySource(hUold);
    }
}

/*
void DiracGraphene1D::ForwardTimeOperator() {
    for( int i = 1; i <= Nx -2; i++ ){

        float flx_x_old, flx_y_old,hflx_x_old, hflx_y_old, tmp_old;

            flx_x_old = Umain[i].p();
            Umain[i].p() = flx_x_old + Umain[i].d2vx();

            hflx_x_old = HoleUmain[i].p();
            HoleUmain[i].p() = hflx_x_old + HoleUmain[i].d2vx();

            tmp_old = Umain[i].tmp();
            Umain[i].tmp() = tmp_old + Umain[i].d2tmp();



    }
}
*/

void DiracGraphene1D::WriteFluidFile(float t){
    int pos_end = Nx - 1 ;
    int pos_ini = 0;
    if (!isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].v()) ||
        !isfinite(Umain[pos_ini].v())) {
        cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
        exit(EXIT_FAILURE);
    }
//data_preview << t << "\t" << Umain[pos_end/2] << "\n";
    data_preview << t << "\t" << Umain[pos_ini] << "\t" << Umain[pos_end] << "\t"<< HoleUmain[pos_ini] <<"\t"<< HoleUmain[pos_end] << "\n";
}


void DiracGraphene1D::CopyFields() {
    float mass, hmass;
    for (int i = 0; i < Nx; ++i) {
        Den[i]=Umain[i].n();
        mass = DensityToMass(Den[i]);
        Vel[i]=Umain[i].v();
        Tmp[i] =Umain[i].tmp();

        HDen[i]=HoleUmain[i].n();
        hmass = DensityToMass(HDen[i]);
        HVel[i]=HoleUmain[i].v();
    }
}

void DiracGraphene1D::SaveSnapShot() {

    Fluid1D::SaveSnapShot();

    hsize_t dim_atr[1] = { 1 };
    DataSpace atr_dataspace = DataSpace (1, dim_atr );

    int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
    snapshot_step = points_per_period / snapshot_per_period;

    this->CopyFields();

    string str_time = to_string(TimeStepCounter / snapshot_step);
    str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
    string name_dataset = "snapshot_" + str_time;

    float currenttime=static_cast<float>(TimeStepCounter) * dt;

    DataSet dataset_hden = GrpHDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceHDen);
    Attribute atr_step_hden = dataset_hden.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_hden = dataset_hden.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_hden.write(HDen, HDF5FLOAT);
    dataset_hden.close();
    atr_step_hden.write(HDF5INT, &TimeStepCounter);
    atr_time_hden.write(HDF5FLOAT , &currenttime);
    atr_step_hden.close();
    atr_time_hden.close();


    DataSet dataset_hvel_x = GrpHVel->createDataSet(name_dataset, HDF5FLOAT, *DataspaceHVel);
    Attribute atr_step_hvel_x = dataset_hvel_x.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_hvel_x = dataset_hvel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_hvel_x.write(Vel, HDF5FLOAT);
    dataset_hvel_x.close();
    atr_step_hvel_x.write(HDF5INT, &TimeStepCounter);
    atr_time_hvel_x.write(HDF5FLOAT , &currenttime);
    atr_step_hvel_x.close();
    atr_time_hvel_x.close();

}


void DiracGraphene1D::CreateHdf5File() {
    TethysBase::CreateHdf5File();
    GrpHDen = new Group(Hdf5File->createGroup("/Data/HoleDensity" ));
    GrpHVel = new Group(Hdf5File->createGroup("/Data/HoleVelocity" ));
}


float DiracGraphene1D::DensityToMass(float density) {
    return sqrt(density*density*density);
}


//____________________ Electrons Functions _________________________


float DiracGraphene1D::EleDensitySource(StateVec1D Uelec, StateVec1D Uholes) {
    //return A*(1.0f - Uelec.n()) + B*(Uelec.n()*Uholes.n() - 1.0f);
    return 0.0f;
}

float DiracGraphene1D::EleVelocitySource(StateVec1D Uelec) {
    //return -1.0f * col_freq * Uelec.v();
    return 0.0f;
}

float DiracGraphene1D::EleDensityFlux(StateVec1D Uelec) {
   return DensityFlux(Uelec);
}

float DiracGraphene1D::EleVelocityFlux(StateVec1D Uelec) {
    return 0.25f * Uelec.v() * Uelec.v() + vel_fer * vel_fer * 0.5f * log(Uelec.n()+1E-6f) + 2.0f * Uelec.S() * Uelec.S() * sqrt(Uelec.n()) - kin_vis*Uelec.grad_v();
}


//____________________ Holes Functions _________________________


float DiracGraphene1D::HolDensitySource(StateVec1D Uelec , StateVec1D Uholes) {
    //return A*(1.0f - Uholes.n()) + B*(Uelec.n()*Uholes.n() - 1.0f);
    return 0.0f;
}

float DiracGraphene1D::HolVelocitySource(StateVec1D Uholes) {
    //return -1.0f * col_freq * Uholes.v();
    return 0.0f;
}

float DiracGraphene1D::HolDensityFlux(StateVec1D Uholes) {
    return DensityFlux(Uholes);
}

float DiracGraphene1D::HolVelocityFlux(StateVec1D Uholes) {
    return 0.25f * Uholes.v() * Uholes.v() + vel_fer * vel_fer * 0.5f * log(Uholes.n()+1E-6f) + 2.0f * Uholes.S() * Uholes.S() * sqrt(Uholes.n()) - kin_vis*Uholes.grad_v();
}






DiracGraphene1D::~DiracGraphene1D(){
    delete[] Den;
    delete[] Vel;
    delete[] Cur;
    delete[] HDen;
    delete[] HVel;
    delete[] vel_snd_arr;
    delete[] HoleUmain;
    delete[] HoleUmid;
}

void DiracGraphene1D::ComputeElectricPotencial() {
//TODO implementar a funçao nao local para calculo do potencial electrico em cada ponto
}

