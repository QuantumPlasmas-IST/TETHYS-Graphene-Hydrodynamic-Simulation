//
// Created by pedro on 05-05-2023.
//

#include "includes/TripleGraphene1DLib.h"


TripleGraphene1D::TripleGraphene1D(SetUpParametersCNP &input_parameters) : DiracGraphene1D(input_parameters) {

    hsize_t dimsf[1];
    dimsf[0] = static_cast<hsize_t>(Nx);
    DataspaceBDen = new DataSpace(RANK, dimsf );
    DataspaceBVel = new DataSpace(RANK, dimsf );
    GrpBDen = nullptr;
    GrpBVel = nullptr;

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
    BDen = new float[Nx];
    BVel = new float[Nx];

    BUmain = new StateVec1D[Nx]();
    BUmid = new StateVec1D[Nx-1]();
}


void TripleGraphene1D::SetSound(){
    for(int i = 0; i<Nx  ;i++){
        vel_snd_arr[i]= vel_snd;
        Umain[i].S()=vel_snd;
        HoleUmain[i].S()=vel_snd;
        BUmain[i].S()=vel_snd;
    }
    for(int i = 0; i<Nx-1  ;i++){
        Umid[i].S()=vel_snd;
        HoleUmid[i].S()=vel_snd;
        BUmid[i].S()=vel_snd;
    }
}


void TripleGraphene1D::InitialCondRand(){
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

        noise =  (float) rd()/maxrand ; //(float) rand()/ (float) RAND_MAX ;
        //HDen[c] = 1.0f + 0.05f * (noise - 0.5f);
        BUmain[c].n()=1.0f + 0.05f * (noise - 0.5f);
    }
    this -> SetSound();
}


void TripleGraphene1D::RichtmyerStep1() {

    DiracGraphene1D::ComputeElectricPotencial("MainGrid");

    for( int i = 0; i <= Nx -2; i++){

        //TODO VERIFICAR EQUACOES E IMPLEMENTAR PARA BOSON

        StateVec1D eUavg{};
        eUavg = 0.5f*(Umain[i+1] + Umain[i]);

        StateVec1D hUavg{};
        hUavg = 0.5f*(HoleUmain[i+1] + HoleUmain[i]);

        StateVec1D bUavg{};
        bUavg = 0.5f*(BUmain[i+1] + BUmain[i]);

        //____________________ Electrons _________________________

        Umid[i].n() = eUavg.n() - 0.5f*(dt/dx)*(EleDensityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - EleDensityFlux(BUmain[i], Umain[i], HoleUmain[i]))
                      + (0.5f*dt) * EleDensitySource(bUavg, eUavg, hUavg) ;

        Umid[i].v() = eUavg.v() - 0.5f*(dt/dx)*(EleVelocityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - EleVelocityFlux(BUmain[i], Umain[i], HoleUmain[i]))
                      + (0.5f*dt) * EleVelocitySource(bUavg, eUavg, hUavg) ;

        //____________________ Holes _____________________________

        HoleUmid[i].n() = hUavg.n() - 0.5f*(dt/dx)*(HolDensityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - HolDensityFlux(BUmain[i], Umain[i], HoleUmain[i]))
                          + (0.5f*dt) * HolDensitySource(bUavg, eUavg, hUavg) ;

        HoleUmid[i].v() = hUavg.v() - 0.5f*(dt/dx)*(HolVelocityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - HolVelocityFlux(BUmain[i], Umain[i], HoleUmain[i]))
                          + (0.5f*dt) * HolVelocitySource(bUavg, eUavg, hUavg) ;

        //____________________ Boson _____________________________

        BUmid[i].n() = bUavg.n() - 0.5f*(dt/dx)*(BDensityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - BDensityFlux(BUmain[i], Umain[i], HoleUmain[i]))
                       + (0.5f*dt) * BDensitySource(bUavg, eUavg, hUavg) ;

        BUmid[i].v() = bUavg.v() - 0.5f*(dt/dx)*(BVelocityFlux(BUmain[i+1], Umain[i+1], HoleUmain[i+1]) - BVelocityFlux(BUmain[i], Umain[i], BUmain[i]))
                       + (0.5f*dt) * BVelocitySource(bUavg, eUavg, hUavg) ;
    }
}

void TripleGraphene1D::RichtmyerStep2() {

    DiracGraphene1D::ComputeElectricPotencial("MidGrid");

    for( int i = 1; i <= Nx -2; i++ ){

        //TODO VERIFICAR EQUACOES E IMPLEMENTAR PARA BOSON

        StateVec1D eUold(Umain[i]);
        StateVec1D hUold(HoleUmain[i]);
        StateVec1D bUold(BUmain[i]);

        //____________________ Electrons _________________________

        Umain[i].n() = eUold.n() - (dt/dx) * (EleDensityFlux(BUmid[i], Umid[i], HoleUmid[i]) - EleDensityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                       + dt * EleDensitySource(bUold, eUold,hUold);

        Umain[i].v() = eUold.v() - (dt/dx)*(EleVelocityFlux(BUmid[i], Umid[i], HoleUmid[i]) - EleVelocityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                       + dt * EleVelocitySource(bUold, eUold,hUold);

        //____________________ Holes _____________________________

        HoleUmain[i].n() = hUold.n() - (dt/dx) * (HolDensityFlux(BUmid[i], Umid[i], HoleUmid[i]) - HolDensityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                           + dt * HolDensitySource(bUold, eUold,hUold);

        HoleUmain[i].v() = hUold.v() - (dt/dx)*(HolVelocityFlux(BUmid[i], Umid[i], HoleUmid[i]) - HolVelocityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                           + dt * HolVelocitySource(bUold, eUold,hUold);

        //____________________ Holes _____________________________

        BUmain[i].n() = bUold.n() - (dt/dx) * (BDensityFlux(BUmid[i], Umid[i], HoleUmid[i]) - BDensityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                        + dt * BDensitySource(bUold, eUold,hUold);

        BUmain[i].v() = bUold.v() - (dt/dx)*(BVelocityFlux(BUmid[i], Umid[i], HoleUmid[i]) - BVelocityFlux(BUmid[i-1], Umid[i-1], HoleUmid[i-1]))
                        + dt * BVelocitySource(bUold, eUold,hUold);

    }
}

void TripleGraphene1D::WriteFluidFile(float t){
    int pos_end = Nx - 1 ;
    int pos_ini = 0;
    if (!isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].v()) ||
        !isfinite(Umain[pos_ini].v())) {
        cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
        exit(EXIT_FAILURE);
    }
//data_preview << t << "\t" << Umain[pos_end/2] << "\n";
    data_preview << t << "\t" << Umain[pos_ini] << "\t" << Umain[pos_end] << "\t"<< HoleUmain[pos_ini] <<"\t"<< HoleUmain[pos_end] << "\t"<< BUmain[pos_ini] <<"\t"<< BUmain[pos_end] <<"\n";
}

void TripleGraphene1D::CopyFields() {
    //float mass, hmass, bmass;
    for (int i = 0; i < Nx; ++i) {
        Den[i]=Umain[i].n();
        //mass = DensityToMass(Den[i]);
        Vel[i]=Umain[i].v();
        //Tmp[i] =Umain[i].tmp();

        HDen[i]=HoleUmain[i].n();
        //hmass = DensityToMass(HDen[i]);
        HVel[i]=HoleUmain[i].v();

        BDen[i] = BUmain[i].n();
        //bmass = DensityToMass(BDen[i]);
        BVel[i] = BUmain[i].v();
    }
}

void TripleGraphene1D::SaveSnapShot() {

    DiracGraphene1D::SaveSnapShot();

    hsize_t dim_atr[1] = { 1 };
    DataSpace atr_dataspace = DataSpace (1, dim_atr );

    int points_per_period = static_cast<int>((2.0 * MAT_PI / this->RealFreq()) / dt);
    snapshot_step = points_per_period / snapshot_per_period;

    this->CopyFields();

    string str_time = to_string(TimeStepCounter / snapshot_step);
    str_time.insert(str_time.begin(), 5 - str_time.length(), '0');
    string name_dataset = "snapshot_" + str_time;

    float currenttime=static_cast<float>(TimeStepCounter) * dt;

    DataSet dataset_bden = GrpBDen->createDataSet(name_dataset, HDF5FLOAT, *DataspaceBDen);
    Attribute atr_step_bden = dataset_bden.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_bden = dataset_bden.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_bden.write(BDen, HDF5FLOAT);
    dataset_bden.close();
    atr_step_bden.write(HDF5INT, &TimeStepCounter);
    atr_time_bden.write(HDF5FLOAT , &currenttime);
    atr_step_bden.close();
    atr_time_bden.close();


    DataSet dataset_bvel_x = GrpBVel->createDataSet(name_dataset, HDF5FLOAT, *DataspaceBVel);
    Attribute atr_step_bvel_x = dataset_bvel_x.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_bvel_x = dataset_bvel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_bvel_x.write(BVel, HDF5FLOAT);
    dataset_bvel_x.close();
    atr_step_bvel_x.write(HDF5INT, &TimeStepCounter);
    atr_time_bvel_x.write(HDF5FLOAT , &currenttime);
    atr_step_bvel_x.close();
    atr_time_bvel_x.close();

}

void TripleGraphene1D::CreateHdf5File() {
    DiracGraphene1D::CreateHdf5File();
    GrpBDen = new Group(Hdf5File->createGroup("/Data/BosonDensity"));
    GrpBVel = new Group(Hdf5File->createGroup("/Data/BosonVelocity"));
}



//____________________ Electron Functions _________________________


float TripleGraphene1D::EleDensitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return -A*(Uboson.n()-Uelec.n());
    //return 0.0f;
}

float TripleGraphene1D::EleVelocitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return (A-B)*DensityFlux(Uelec)+12.0f*A*DensityFlux(Uboson);
    //return 0.0f;
}

float TripleGraphene1D::EleDensityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return DensityFlux(Uelec);
    //return 0.0f;
}

float TripleGraphene1D::EleVelocityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return 0.25f * Uelec.v() * Uelec.v() + vel_fer * vel_fer * Uelec.n()+ 2.0f * Uelec.S() * Uelec.S() * Uelec.phi() + vel_pm * vel_pm * Uboson.n();// - kin_vis*Uelec.grad_v();
    //return 0.0f;
}


//____________________ Hole Functions _________________________


float TripleGraphene1D::HolDensitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return -A*(Uboson.n()-Uholes.n());
    return 0.0f;
}

float TripleGraphene1D::HolVelocitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return (A-B)*DensityFlux(Uholes)+12.0f*A*DensityFlux(Uboson);
    //return 0.0f;
}

float TripleGraphene1D::HolDensityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    //return 0.0f;
    return DensityFlux(Uholes);
}

float TripleGraphene1D::HolVelocityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    //return 0.0f;
    return 0.25f * Uholes.v() * Uholes.v() + vel_fer * vel_fer * Uholes.n() - 2.0f * Uholes.S() * Uholes.S() * Uholes.phi() + vel_pm *vel_pm * Uboson.n();// - kin_vis*Uholes.grad_v();
}

//____________________ Boson Functions _________________________

float TripleGraphene1D::BDensitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return -A*(Uelec.n()+Uholes.n()-2.0f*Uboson.n());
    //return 0.0f;
}

float TripleGraphene1D::BVelocitySource(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    return (1/12.0f)*A*(DensityFlux(Uelec)+ DensityFlux(Uholes))+2.0f*A* DensityFlux(Uboson);
    //return 0.0f;
}

float TripleGraphene1D::BDensityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    //return 0.0f;
    return DensityFlux(Uboson);
}

float TripleGraphene1D::BVelocityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes) {
    //return 0.0f;
    return 0.25f * Uboson.v() * Uboson.v() - vel_p * vel_p * (Uholes.n() + Uelec.n());
}


void TripleGraphene1D::InitialCondTest(){
    for (int i = 0; i < Nx; i++ ){

        //Umain[i].v()= -1.5f; // (i>3*Nx/8 && i<5*Nx/8 ) ? 3.0f : 0.0f; //1.5f;//
        //HoleUmain[i].v()= 1.5f;//1.0f/(1.0f+5.0f* pow(cosh((i*dx-0.5f)*12.0f),2.f));
        HoleUmain[i].n()= 0.2f + 0.1f/ pow(cosh((i*dx-0.5f)*12.0f),2.f); //(i>3*Nx/8 && i<5*Nx/8 ) ? 1.0f : 0.1f; //0.2f+0.2f/ pow(cosh((i*dx-0.5f)*12.0f),2);//
        Umain[i].n()= 0.2f + 0.1f/ pow(cosh((i*dx-0.5f)*12.0f),2.f); //(i>3*Nx/8 && i<5*Nx/8 ) ? 1.0f : 0.1f; //0.2f+0.2f/ pow(cosh((i*dx-0.5f)*12.0f),2);//
        BUmain[i].n() = 0.2f;//+0.1f/ pow(cosh((i*dx-0.5f)*12.0f),2.f);

        //HoleUmain[i].n()=0.05f* cos((2*MAT_PI/0.2f)*i*dx)+0.25f;
        //Umain[i].n()=0.05f* cos((2*MAT_PI/0.2f)*i*dx)+0.25f;//
    }
    this->SetSound();
}

TripleGraphene1D::~TripleGraphene1D(){
    delete[] BDen;
    delete[] BVel;
    delete[] BUmain;
    delete[] BUmid;
}