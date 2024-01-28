//
// Created by pedro on 07-01-2024.
//

#include "TripleGraphene2DLib.h"


TripleGraphene2D::TripleGraphene2D(SetUpParametersCNP &input_parameters) : DiracGraphene2D(input_parameters) {

    hsize_t dimsf[2];
    dimsf[0] = static_cast<hsize_t>(Ny);
    dimsf[1] = static_cast<hsize_t>(Nx);
    DataspaceBDen = new DataSpace(RANK, dimsf );
    DataspaceBVelX = new DataSpace(RANK, dimsf );
    DataspaceBVelY = new DataSpace(RANK, dimsf );
    GrpBDen = nullptr;
    GrpBVelX = nullptr;
    GrpBVelY = nullptr;

    vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
    col_freq = input_parameters.CollisionFrequency ; // collision_frequency
    cyc_freq = input_parameters.CyclotronFrequency ; //cyclotron_frequency
    therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity

    vel_therm = input_parameters.ThermalVelocity ;
    A = input_parameters.Diffusive_sourceterm ;
    B = input_parameters.Creation_sourceterm ;

    char buffer [100];
    sprintf (buffer, "S=%.2fvF=%.2fvT=%.2fA=%.3fB=%.3fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, vel_therm, A, B, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
    file_infix = buffer;

    // main grid variables Nx*Ny
    BDen 		= new float[Nx * Ny]();
    BVelX 		= new float[Nx * Ny]();
    BVelY 		= new float[Nx * Ny]();

    BUmain = new StateVec2D[Nx * Ny]();
    BUmid = new StateVec2D[(Nx - 1) * (Ny - 1)]();
}

void TripleGraphene2D::SetSound(){
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

float TripleGraphene2D::FermionDensityToMass(float density) {
    return sqrt(density);
}

float TripleGraphene2D::BosonDensityToMass(float density) {
    return pow(density, float (0.75));
}

void TripleGraphene2D::InitialCondRand(){
    random_device rd;
    float maxrand;
    maxrand = (float) random_device::max();
    for (int c = 0; c < Nx*Ny; c++ ){
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

void TripleGraphene2D::InitialCondPointDen() {
    for (int c = 0; c < Nx*Ny; c++ ){
        //Den[c] = 1.0f;
        //HDen[c] = 1.0f;
        Umain[c].n()=1.0f;
        HoleUmain[c].n()=1.0f;
        BUmain[c].n()=1.0f;
    }
    for (int c = 3*Ny/7; c < 4*Ny/7; c++){
        for (int g = 3*Nx/7; g < 4*Nx/7; g++){
            Umain[c*Nx+g].n()=1.2f;
            HoleUmain[c*Nx+g].n()=1.2f;
            //BUmain[c*Nx+g].n()=1.02f;
        }
    }
    for (int c = 3*Ny/7; c < 4*Ny/7; c++){
        for (int g = 4*Nx/7; g < 5*Nx/7; g++){
            //Den[c*Nx+g] = 0.8f;
            //HDen[c*Nx+g] = 1.2f;
            //Umain[c*Nx+g].n()=0.8f;
            //HoleUmain[c*Nx+g].n()=1.2f;
            BUmain[c*Nx+g].n()=1.1f;
        }
    }
    this->SetSound();
}



void TripleGraphene2D::InitialCondTest() {

    float periods = 1.0f;

    for (int g = 0; g < Nx ; g++) { //g coordenada x
        for (int c = 0; c < Ny ; c++){ //c coordenada y

            /*
            BUmain[c*Nx+g].n()= 0.5f + 0.5f / pow(cosh(((g-Nx/4.f) * (g-Nx/4.f) + (c-Ny/4.f) * (c-Ny/4.f)) * 0.06f *0.06f),2.f);// + 0.5f / pow(cosh(((g-3.0f*Nx/4.f) * (g-3.0f*Nx/4.f) + (c-3.0f*Ny/4.f) * (c-3.0f*Ny/4.f)) * 0.06f *0.06f),2.f);
            HoleUmain[c*Nx+g].n() = 1.f;
            //HoleUmain[c*Nx+g].n()= 1.f + 0.1f / pow(cosh(((g-Nx/2) * (g-Nx/2) + (c-Ny/2) * (c-Ny/2)) * .04f*0.04f),2.f);
            Umain[c*Nx+g].n()=1.f;// + 0.1f / pow(cosh(((g-Nx/2) * (g-Nx/2) + (c-Ny/2) * (c-Ny/2)) * 0.1f*0.1f),2.f);
            */

            /*
            if ( pow(g - Nx/4.f,2.0f)+pow(c - Ny/4.f,2.f) <= 25.f*25.f){
                BUmain[c*Nx+g].px() = 10.f;// * BosonDensityToMass(BUmain[c*Nx+g].n());// * BUmain[c*Nx+g].n();
            }
            */

            /*
            if ( pow(g - 3.f*Nx/4.f  ,2.0f)+pow(c - 3*Ny/4.f,2.f) <= 25.f*25.f){
                BUmain[c*Nx+g].px() = - 10.f* BosonDensityToMass(BUmain[c*Nx+g].n());// * BUmain[c*Nx+g].n();
            }
            */


            Umain[c*Nx+g].n()= 0.2f * cos( periods * 2*MAT_PI * (g-Nx/2) / Nx ) * cos(periods * 2*MAT_PI * (c-Ny/2) / Ny ) + 1.0f;  //condição periodica em x e y
            HoleUmain[c*Nx+g].n() = Umain[c*Nx+g].n();
            BUmain[c*Nx+g].n() = 1.f;

        }
    }
    this -> SetSound();
}


void TripleGraphene2D::ChooseGridPointers(const string &grid) {
    if(grid == "MidGrid"){
        ptr_StateVec = Umain;
        ptr_StateVecHole = HoleUmain;
        ptr_StateVecBoson = BUmain;
    }if(grid == "MainGrid"){
        ptr_StateVec = Umid;
        ptr_StateVecHole = HoleUmid;
        ptr_StateVecBoson = BUmid;
    }
}



void TripleGraphene2D::RichtmyerStep1(){

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ChooseGridPointers("MidGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,Den,FlxX,FlxY,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
    for(int ks=0; ks<=Nx*Ny-Nx-Ny; ks++){ //correr todos os pontos da grelha secundaria de den_mid
        GridPoint2D midpoint(ks, Nx, Ny, true);

        //electrons
        StateVec2D eUavg{};
        eUavg = 0.25f * (Umain[midpoint.SW] + Umain[midpoint.SE] + Umain[midpoint.NW] + Umain[midpoint.NE]);
        StateVec2D eUNorth{};
        StateVec2D eUSouth{};
        StateVec2D eUEast{};
        StateVec2D eUWest{};
        eUNorth = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.NW]);
        eUSouth = 0.5f*(Umain[midpoint.SE]+Umain[midpoint.SW]);
        eUEast = 0.5f*(Umain[midpoint.NE]+Umain[midpoint.SE]);
        eUWest = 0.5f*(Umain[midpoint.NW]+Umain[midpoint.SW]);
        //holes
        StateVec2D hUavg{};
        hUavg = 0.25f * (HoleUmain[midpoint.SW] + HoleUmain[midpoint.SE] + HoleUmain[midpoint.NW] + HoleUmain[midpoint.NE]);
        StateVec2D hUNorth{};
        StateVec2D hUSouth{};
        StateVec2D hUEast{};
        StateVec2D hUWest{};
        hUNorth = 0.5f*(HoleUmain[midpoint.NE]+HoleUmain[midpoint.NW]);
        hUSouth = 0.5f*(HoleUmain[midpoint.SE]+HoleUmain[midpoint.SW]);
        hUEast = 0.5f*(HoleUmain[midpoint.NE]+HoleUmain[midpoint.SE]);
        hUWest = 0.5f*(HoleUmain[midpoint.NW]+HoleUmain[midpoint.SW]);
        //bosons
        StateVec2D bUavg{};
        bUavg = 0.25f * (BUmain[midpoint.SW] + BUmain[midpoint.SE] + BUmain[midpoint.NW] + BUmain[midpoint.NE]);
        StateVec2D bUNorth{};
        StateVec2D bUSouth{};
        StateVec2D bUEast{};
        StateVec2D bUWest{};
        bUNorth = 0.5f*(BUmain[midpoint.NE]+BUmain[midpoint.NW]);
        bUSouth = 0.5f*(BUmain[midpoint.SE]+BUmain[midpoint.SW]);
        bUEast = 0.5f*(BUmain[midpoint.NE]+BUmain[midpoint.SE]);
        bUWest = 0.5f*(BUmain[midpoint.NW]+BUmain[midpoint.SW]);

        //____________________ Electrons _________________________

        Umid[ks].n() =  eUavg.n()
                        -0.5f*(dt/dx)*(EleDensityFluxX(eUEast,hUEast,bUEast) - EleDensityFluxX(eUWest,hUWest,bUWest))
                        -0.5f*(dt/dy)*(EleDensityFluxY(eUNorth,hUNorth,bUNorth) - EleDensityFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt* EleDensitySource(eUavg,hUavg,bUavg);

        Umid[ks].px() = eUavg.px()
                        -0.5f*(dt/dx)*(EleXMomentumFluxX(eUEast,hUEast,bUEast) - EleXMomentumFluxX(eUWest,hUWest,bUWest))
                        -0.5f*(dt/dy)*(EleXMomentumFluxY(eUNorth,hUNorth,bUNorth) - EleXMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*EleXMomentumSource(eUavg,hUavg,bUavg);

        Umid[ks].py() = eUavg.py()
                        -0.5f*(dt/dx)*(EleYMomentumFluxX(eUEast,hUEast,bUEast) - EleYMomentumFluxX(eUWest,hUWest,bUWest))
                        -0.5f*(dt/dy)*(EleYMomentumFluxY(eUNorth,hUNorth,bUNorth) - EleYMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*EleYMomentumSource(eUavg,hUavg,bUavg);

        //________________________________ HOLE ____________________________

        HoleUmid[ks].n() =  hUavg.n()
                            -0.5f*(dt/dx)*(HolDensityFluxX(eUEast,hUEast,bUEast) - HolDensityFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(HolDensityFluxY(eUNorth,hUNorth,bUNorth) - HolDensityFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt* HolDensitySource(eUavg,hUavg,bUavg);

        HoleUmid[ks].px() = hUavg.px()
                            -0.5f*(dt/dx)*(HolXMomentumFluxX(eUEast,hUEast,bUEast) - HolXMomentumFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(HolXMomentumFluxY(eUNorth,hUNorth,bUNorth) - HolXMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*HolXMomentumSource(eUavg,hUavg,bUavg);

        HoleUmid[ks].py() = hUavg.py()
                            -0.5f*(dt/dx)*(HolYMomentumFluxX(eUEast,hUEast,bUEast) - HolYMomentumFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(HolYMomentumFluxY(eUNorth,hUNorth,bUNorth) - HolYMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*HolYMomentumSource(eUavg,hUavg,bUavg);
        //________________________________ BOSON ____________________________

        BUmid[ks].n() =  bUavg.n()
                            -0.5f*(dt/dx)*(BosonDensityFluxX(eUEast,hUEast,bUEast) - BosonDensityFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(BosonDensityFluxY(eUNorth,hUNorth,bUNorth) - BosonDensityFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt* BosonDensitySource(eUavg,hUavg,bUavg);

        BUmid[ks].px() = bUavg.px()
                            -0.5f*(dt/dx)*(BosonXMomentumFluxX(eUEast,hUEast,bUEast) - BosonXMomentumFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(BosonXMomentumFluxY(eUNorth,hUNorth,bUNorth) - BosonXMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*BosonXMomentumSource(eUavg,hUavg,bUavg);

        BUmid[ks].py() = bUavg.py()
                            -0.5f*(dt/dx)*(BosonYMomentumFluxX(eUEast,hUEast,bUEast) -BosonYMomentumFluxX(eUWest,hUWest,bUWest))
                            -0.5f*(dt/dy)*(BosonYMomentumFluxY(eUNorth,hUNorth,bUNorth) - BosonYMomentumFluxY(eUSouth,hUSouth,bUSouth))+0.5f*dt*BosonYMomentumSource(eUavg,hUavg,bUavg);
    }
}

void TripleGraphene2D::RichtmyerStep2(){

    ChooseGridPointers("MainGrid");
//#pragma omp parallel for default(none) shared(Nx,Ny,dt,dx,dy,FlxX,FlxY,Den,Tmp,den_dx,den_dy,ptr_den,ptr_px,ptr_py,ptr_snd,ptr_tmp,ptr_velXdx,ptr_velXdy,ptr_velYdx,ptr_velYdy,ptr_dendx,ptr_dendy)
    for(int kp=1+Nx; kp<=Nx*Ny-Nx-2; kp++){ //correr a grelha principal evitando as fronteiras
        GridPoint2D mainpoint(kp, Nx, Ny, false);
        if( kp%Nx!=Nx-1 && kp%Nx!=0){

            //electrons
            StateVec2D eUold(Umain[kp]);
            StateVec2D eUNorth{};
            StateVec2D eUSouth{};
            StateVec2D eUEast{};
            StateVec2D eUWest{};
            eUNorth = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.NW]);
            eUSouth = 0.5f*(Umid[mainpoint.SE]+Umid[mainpoint.SW]);
            eUEast = 0.5f*(Umid[mainpoint.NE]+Umid[mainpoint.SE]);
            eUWest = 0.5f*(Umid[mainpoint.NW]+Umid[mainpoint.SW]);

            //holes
            StateVec2D hUold(HoleUmain[kp]);
            StateVec2D hUNorth{};
            StateVec2D hUSouth{};
            StateVec2D hUEast{};
            StateVec2D hUWest{};
            hUNorth = 0.5f*(HoleUmid[mainpoint.NE]+HoleUmid[mainpoint.NW]);
            hUSouth = 0.5f*(HoleUmid[mainpoint.SE]+HoleUmid[mainpoint.SW]);
            hUEast = 0.5f*(HoleUmid[mainpoint.NE]+HoleUmid[mainpoint.SE]);
            hUWest = 0.5f*(HoleUmid[mainpoint.NW]+HoleUmid[mainpoint.SW]);

            //boson
            StateVec2D bUold(BUmain[kp]);
            StateVec2D bUNorth{};
            StateVec2D bUSouth{};
            StateVec2D bUEast{};
            StateVec2D bUWest{};
            bUNorth = 0.5f*(BUmid[mainpoint.NE]+BUmid[mainpoint.NW]);
            bUSouth = 0.5f*(BUmid[mainpoint.SE]+BUmid[mainpoint.SW]);
            bUEast = 0.5f*(BUmid[mainpoint.NE]+BUmid[mainpoint.SE]);
            bUWest = 0.5f*(BUmid[mainpoint.NW]+BUmid[mainpoint.SW]);


            //____________________ Electrons _________________________


            Umain[kp].n() = eUold.n()
                            - (dt/dx)*(EleDensityFluxX(eUEast,hUEast,bUEast) - EleDensityFluxX(eUWest,hUWest,bUWest))
                            - (dt/dy)*(EleDensityFluxY(eUNorth,hUNorth,bUNorth) - EleDensityFluxY(eUSouth,hUSouth,bUSouth))+ dt*EleDensitySource(eUold,hUold,bUold);
            Umain[kp].px() = eUold.px()
                             - (dt/dx)*(EleXMomentumFluxX(eUEast,hUEast,bUEast) - EleXMomentumFluxX(eUWest,hUWest,bUWest))
                             - (dt/dy)*(EleXMomentumFluxY(eUNorth,hUNorth,bUNorth) - EleXMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*EleXMomentumSource(eUold,hUold,bUold);

            Umain[kp].py() = eUold.py()
                             - (dt/dx)*(EleYMomentumFluxX(eUEast,hUEast,bUEast) - EleYMomentumFluxX(eUWest,hUWest,bUWest))
                             - (dt/dy)*(EleYMomentumFluxY(eUNorth,hUNorth,bUNorth) - EleYMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*EleYMomentumSource(eUold,hUold,bUold);


            //________________________________ HOLE ____________________________

            HoleUmain[kp].n() = hUold.n()
                                - (dt/dx)*(HolDensityFluxX(eUEast,hUEast,bUEast) - HolDensityFluxX(eUWest,hUWest,bUWest))
                                - (dt/dy)*(HolDensityFluxY(eUNorth,hUNorth,bUNorth) - HolDensityFluxY(eUSouth,hUSouth,bUSouth))+ dt*HolDensitySource(eUold,hUold,bUold);
            HoleUmain[kp].px() = hUold.px()
                                 - (dt/dx)*(HolXMomentumFluxX(eUEast,hUEast,bUEast) - HolXMomentumFluxX(eUWest,hUWest,bUWest))
                                 - (dt/dy)*(HolXMomentumFluxY(eUNorth,hUNorth,bUNorth) - HolXMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*HolXMomentumSource(eUold,hUold,bUold);

            HoleUmain[kp].py() = hUold.py()
                                 - (dt/dx)*(HolYMomentumFluxX(eUEast,hUEast,bUEast) - HolYMomentumFluxX(eUWest,hUWest,bUWest))
                                 - (dt/dy)*(HolYMomentumFluxY(eUNorth,hUNorth,bUNorth) - HolYMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*HolYMomentumSource(eUold,hUold,bUold);

            //________________________________ BOSON ____________________________

            BUmain[kp].n() = bUold.n()
                                - (dt/dx)*(BosonDensityFluxX(eUEast,hUEast,bUEast) - BosonDensityFluxX(eUWest,hUWest,bUWest))
                                - (dt/dy)*(BosonDensityFluxY(eUNorth,hUNorth,bUNorth) - BosonDensityFluxY(eUSouth,hUSouth,bUSouth))+ dt*BosonDensitySource(eUold,hUold,bUold);
            BUmain[kp].px() = bUold.px()
                                 - (dt/dx)*(BosonXMomentumFluxX(eUEast,hUEast,bUEast) - BosonXMomentumFluxX(eUWest,hUWest,bUWest))
                                 - (dt/dy)*(BosonXMomentumFluxY(eUNorth,hUNorth,bUNorth) - BosonXMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*BosonXMomentumSource(eUold,hUold,bUold);

            BUmain[kp].py() = bUold.py()
                                 - (dt/dx)*(BosonYMomentumFluxX(eUEast,hUEast,bUEast) - BosonYMomentumFluxX(eUWest,hUWest,bUWest))
                                 - (dt/dy)*(BosonYMomentumFluxY(eUNorth,hUNorth,bUNorth) - BosonYMomentumFluxY(eUSouth,hUSouth,bUSouth))+ dt*BosonYMomentumSource(eUold,hUold,bUold);
        }
    }
}


float TripleGraphene2D::EleDensityFluxX(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float den=Uelec.n();
    float mass = FermionDensityToMass(den);
    float px=Uelec.px();
    return px / mass;
}

float TripleGraphene2D::EleDensityFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float den=Uelec.n();
    float mass = FermionDensityToMass(den);
    float py=Uelec.py();
    return py / mass;
}

float TripleGraphene2D::EleXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {

    float sound=vel_snd;
    float den=Uelec.n();
    float vel_pm = 0.001f;
    float hden=Uholes.n();
    float px=Uelec.px();
    float mass=FermionDensityToMass(den);

    return px * px / (den * mass) + sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden) - vel_pm * vel_pm * Uboson.n();
}

float TripleGraphene2D::EleXMomentumFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float den=Uelec.n();
    float px=Uelec.px();
    float py=Uelec.py();
    float mass=FermionDensityToMass(den);
    return px * py / (den * mass);
}


float TripleGraphene2D::EleYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float sound=vel_snd;
    float den=Uelec.n();
    float hden=Uholes.n();
    float vel_pm = 0.001f;
    float py=Uelec.py();
    float mass=FermionDensityToMass(den);
    return py * py / (den * mass) + sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden) - vel_pm * vel_pm * Uboson.n();
}

float TripleGraphene2D::EleYMomentumFluxX(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float den=Uelec.n();
    float px=Uelec.px();
    float py=Uelec.py();
    float mass=FermionDensityToMass(den);

    return px * py / (den * mass);
}


float TripleGraphene2D::EleDensitySource(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::EleXMomentumSource( __attribute__((unused))StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::EleYMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}


//________________________________ HOLE ____________________________

float TripleGraphene2D::HolDensityFluxX(__attribute__((unused))  StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float hden=Uholes.n();
    float px=Uholes.px();
    float mass = FermionDensityToMass(hden);
    return px / mass;
}

float TripleGraphene2D::HolDensityFluxY(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float hden=Uholes.n();
    float py=Uholes.py();
    float mass = FermionDensityToMass(hden);
    return py / mass;
}

float TripleGraphene2D::HolXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float vel_pm = 0.001f;
    float sound=vel_snd;
    float den=Uelec.n();
    float hden=Uholes.n();
    float px=Uholes.px();
    float mass=FermionDensityToMass(hden);

    return px * px / (hden * mass) - sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden) - vel_pm * vel_pm * Uboson.n();
}

float TripleGraphene2D::HolXMomentumFluxY(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float hden=Uholes.n();
    float px=Uholes.px();
    float py=Uholes.py();
    float mass=FermionDensityToMass(hden);
    return px * py / (hden * mass);
}


float TripleGraphene2D::HolYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float sound=vel_snd;
    float vel_pm = 0.001f;
    float den=Uelec.n();
    float hden=Uholes.n();
    float py=Uholes.py();
    float mass=FermionDensityToMass(hden);
    return py * py / (hden * mass) - sound * sound * (den - hden) + vel_therm * vel_therm * (den + hden) - vel_pm * vel_pm * Uboson.n();
}

float TripleGraphene2D::HolYMomentumFluxX(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    float hden=Uholes.n();
    float px=Uholes.px();
    float py=Uholes.py();
    float mass=FermionDensityToMass(hden);
    return px * py / (hden * mass);
}

float TripleGraphene2D::HolDensitySource(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::HolXMomentumSource(__attribute__((unused)) StateVec2D Uelec ,__attribute__((unused))  StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::HolYMomentumSource(__attribute__((unused)) StateVec2D Uelec ,__attribute__((unused))  StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

//________________________________ BOSON ____________________________

float TripleGraphene2D::BosonDensitySource(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::BosonXMomentumSource(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::BosonYMomentumSource(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    return 0.0f;
}

float TripleGraphene2D::BosonDensityFluxX(__attribute__((unused)) StateVec2D Uelec, __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float pbx = Uboson.px();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    return pbx / bmass;
}

float TripleGraphene2D::BosonDensityFluxY(__attribute__((unused)) StateVec2D Uelec,__attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson) {
    float pby = Uboson.px();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    return pby / bmass;
}

float TripleGraphene2D::BosonXMomentumFluxX(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    float vel_ref = .001f;
    float pbx = Uboson.px();
    float den = Uelec.n();
    float hden = Uholes.n();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    return (pbx * pbx)/(bden * bmass) - vel_ref * vel_ref * (den + hden);
}

float TripleGraphene2D::BosonXMomentumFluxY(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    float pbx = Uboson.px();
    float pby = Uboson.py();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    return (pbx * pby)/(bden * bmass);
}

float TripleGraphene2D::BosonYMomentumFluxX(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    float pbx = Uboson.px();
    float pby = Uboson.py();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    return (pby * pbx)/(bden * bmass) ;
}

float TripleGraphene2D::BosonYMomentumFluxY(StateVec2D Uelec, StateVec2D Uholes, StateVec2D Uboson) {
    float pby = Uboson.py();
    float bden = Uboson.n();
    float bmass = BosonDensityToMass(bden);
    float vel_ref = .001f;
    float den = Uelec.n();
    float hden = Uholes.n();
    return (pby * pby)/(bden * bmass) - vel_ref * vel_ref * (den + hden);
}

void TripleGraphene2D::WriteFluidFile(float t){
    int j=Ny/2;
    int pos_end = Nx - 1 + j*Nx ;
    int pos_ini = j*Nx ;
    if(!isfinite(Umain[pos_ini].n()) || !isfinite(Umain[pos_end].n()) || !isfinite(Umain[pos_ini].px()) || !isfinite(Umain[pos_end].px())){
        cerr << "ERROR: numerical method failed to converge" <<"\nExiting"<< endl;
        CloseHdf5File();
        exit(EXIT_FAILURE);
    }
    data_preview << t <<"\t"<< Umain[pos_ini] <<"\t"<<Umain[pos_end]<< "\t"<< HoleUmain[pos_ini] <<"\t"<< HoleUmain[pos_end]<< "\t" << BUmain[pos_ini] << "\t" << BUmain[pos_end] << "\n";
}

void TripleGraphene2D::ForwardTimeOperator() {
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
        float flx_x_old, flx_y_old,hflx_x_old, hflx_y_old, bflx_x_old, bflx_y_old,tmp_old;
        if (kp % Nx != Nx - 1 && kp % Nx != 0) {
            flx_x_old = Umain[kp].px();
            flx_y_old = Umain[kp].py();
            Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
            Umain[kp].py() = flx_y_old + Umain[kp].d2vy();

            hflx_x_old = HoleUmain[kp].px();
            hflx_y_old = HoleUmain[kp].py();
            HoleUmain[kp].px() = hflx_x_old + HoleUmain[kp].d2vx();
            HoleUmain[kp].py() = hflx_y_old + HoleUmain[kp].d2vy();

            bflx_x_old = BUmain[kp].px();
            bflx_y_old = BUmain[kp].py();
            BUmain[kp].px() = bflx_x_old + BUmain[kp].d2vx();
            BUmain[kp].py() = bflx_y_old + BUmain[kp].d2vy();

            tmp_old = Umain[kp].tmp();
            Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();


        }
    }
}

void TripleGraphene2D::ForwardTimeOperator(char field) {
    for (int kp = 1 + Nx; kp <= Nx * Ny - Nx - 2; kp++) {
        float flx_x_old, flx_y_old,hflx_x_old, hflx_y_old, bflx_x_old, bflx_y_old, tmp_old;
        if (kp % Nx != Nx - 1 && kp % Nx != 0) {
            switch(field) {
                case 'T':
                    tmp_old = Umain[kp].tmp();
                    Umain[kp].tmp() = tmp_old + Umain[kp].d2tmp();
                    break;
                case 'V':
                    flx_x_old = Umain[kp].px();
                    flx_y_old = Umain[kp].py();
                    Umain[kp].px() = flx_x_old + Umain[kp].d2vx();
                    Umain[kp].py() = flx_y_old + Umain[kp].d2vy();

                    hflx_x_old = HoleUmain[kp].px();
                    hflx_y_old = HoleUmain[kp].py();
                    HoleUmain[kp].px() = hflx_x_old + HoleUmain[kp].d2vx();
                    HoleUmain[kp].py() = hflx_y_old + HoleUmain[kp].d2vy();

                    bflx_x_old = BUmain[kp].px();
                    bflx_y_old = BUmain[kp].py();
                    BUmain[kp].px() = bflx_x_old + BUmain[kp].d2vx();
                    BUmain[kp].py() = bflx_y_old + BUmain[kp].d2vy();

                    break;
                default: ;
            }
        }
    }
}

void TripleGraphene2D::CopyFields() {
    float mass, hmass, bmass;
    for (int i = 0; i < Nx*Ny; ++i) {
        Den[i]=Umain[i].n();
        mass= FermionDensityToMass(Den[i]);
        VelX[i]=Umain[i].px()/mass;
        VelY[i]=Umain[i].py()/mass;
        //Tmp[i] =Umain[i].tmp();

        HDen[i]=HoleUmain[i].n();
        hmass= FermionDensityToMass(HDen[i]);
        HVelX[i]=HoleUmain[i].px()/hmass;
        HVelY[i]=HoleUmain[i].py()/hmass;

        BDen[i]=BUmain[i].n();
        bmass = BosonDensityToMass(BDen[i]);
        BVelX[i]=BUmain[i].px()/bmass;
        BVelY[i]=BUmain[i].py()/bmass;

    }
}

void TripleGraphene2D::SaveSnapShot() {

    DiracGraphene2D::SaveSnapShot();

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


    DataSet dataset_bvel_x = GrpBVelX->createDataSet(name_dataset, HDF5FLOAT, *DataspaceBVelX);
    Attribute atr_step_bvel_x = dataset_bvel_x.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_bvel_x = dataset_bvel_x.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_bvel_x.write(BVelX, HDF5FLOAT);
    dataset_bvel_x.close();
    atr_step_bvel_x.write(HDF5INT, &TimeStepCounter);
    atr_time_bvel_x.write(HDF5FLOAT , &currenttime);
    atr_step_bvel_x.close();
    atr_time_bvel_x.close();

    DataSet dataset_bvel_y = GrpBVelY->createDataSet(name_dataset, HDF5FLOAT, *DataspaceBVelY);
    Attribute atr_step_bvel_y = dataset_bvel_y.createAttribute("time step", HDF5INT, atr_dataspace);
    Attribute atr_time_bvel_y = dataset_bvel_y.createAttribute("time", HDF5FLOAT, atr_dataspace);
    dataset_bvel_y.write(BVelY, HDF5FLOAT);
    dataset_bvel_y.close();
    atr_step_bvel_y.write(HDF5INT, &TimeStepCounter);
    atr_time_bvel_y.write(HDF5FLOAT , &currenttime);
    atr_step_bvel_y.close();
    atr_time_bvel_y.close();

}

void TripleGraphene2D::CreateHdf5File() {
    DiracGraphene2D::CreateHdf5File();
    GrpBDen = new Group(Hdf5File->createGroup("/Data/BosonDensity" ));
    GrpBVelX = new Group(Hdf5File->createGroup("/Data/BosonVelocityX" ));
    GrpBVelY = new Group(Hdf5File->createGroup("/Data/BosonVelocityY"));
}

TripleGraphene2D::~TripleGraphene2D(){
    delete[] BDen;
    delete[] BVelX;
    delete[] BVelY;
    delete[] BUmain;
    delete[] BUmid;
}