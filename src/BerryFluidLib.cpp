/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281                                                                    *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/BerryFluidLib.h"
#include <cmath>

//float theta = .0;
//float A = .01;
float elec_charge=1.60217663*pow(10, -19);
float d0 = 0.001;
float permitivity=8.85418782*pow(10, -12);
float gap = 1.60*pow(10, -20); //100meV
float hbar = 1.05457182*pow(10, -34);
float chem_pot = 0.1;

BerryFluid::BerryFluid(SetUpParameters &input_parameters) : Fluid2D(input_parameters) {
    vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
    col_freq = input_parameters.CollisionFrequency ;
    //therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity
    //theta = input_parameters.ElectricFieldAngle; //angle described by the electric field in relation to xx
    //delta = input_parameters.BandGap; //band gap
    //miu = input_parameters.ChemicalPotential
    char buffer [100];
    sprintf (buffer, "S=%.2fvF=%.2fvis=%.3fodd=%.3fl=%.3f", vel_snd, vel_fer, kin_vis, odd_vis, col_freq);
    file_infix = buffer;
}


float BerryFluid::DensityToMass(float density)
{
    return sqrt(density*density*density);
}

/*
void BerryFluid::SetSimulationTime(){
    float s;
    s=this->GetVelSnd();
    this->SetTmax(5.0f+0.02f*s+20.0f/s);
}*/


/*
float AnomalousStressTensor(float temp, float miu, float delta) {
    float z = exp((miu-1))/(temp/delta);
    float alpha1 = 0.389037*temp/delta;
    float alpha2 = 0.606887*temp/delta;
    float beta1 = 0.620443*temp/delta;
    float beta2 =2.710214*temp/delta;

    float term1 = gap * elec_charge * MAT_PI / (pow(2*MAT_PI*vel_fer, 2) * hbar);
    float term2 = z / (2 * temp/delta);
    float term3 = (temp/delta / (z + 1) - alpha1 /(1 + beta1) * hyp2f1(2, beta1 + 1, beta1 + 2, -z) - alpha2 / (1 + beta2) * hyp2f1(2, beta2 + 1, beta2 + 2, -z));

    return term1 * term2 * term3;
}*/


float BerryFluid::XMomentumSource(StateVec2D U) {
    return -1.0f*col_freq*U.px();
}


float BerryFluid::YMomentumSource(StateVec2D U) {
    return -1.0f*col_freq*U.py();
}


float BerryFluid::XMomentumFluxX(StateVec2D U) {
    float sound=U.S();
    float den=U.n();
    float px=U.px();
    float mass=DensityToMass(den);
    float Vxy=U.dxvy();
    float ny=U.dyn();
    float temp = U.tmp();
    //float D = AnomalousStressTensor(temp, chem_pot, gap);
    float D = 0.1;
    float anomalousXX = px/mass*elec_charge*d0*ny*D/permitivity;

    return px * px / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den + odd_vis*Vxy + anomalousXX;
}

float BerryFluid::XMomentumFluxY(StateVec2D U) {
    float den=U.n();
    float px=U.px();
    float py=U.py();
    float mass=DensityToMass(den);
    float Vyy=U.dyvy();
    float nx=U.dxn();
    float temp = U.tmp();
    //float D = AnomalousStressTensor(temp, chem_pot, gap);
    float D = 0.1;
    float anomalousXY = -px/mass*elec_charge*d0*nx*D/permitivity;

    return px * py / mass + odd_vis*Vyy + anomalousXY;
}

float BerryFluid::YMomentumFluxX(StateVec2D U) {
    float den=U.n();
    float px=U.px();
    float py=U.py();
    float Vxx=U.dxvx();
    float mass=DensityToMass(den);
    float ny=U.dyn();
    float temp = U.tmp();
    //float D = AnomalousStressTensor(temp, chem_pot, gap);
    float D = 0.1;
    float anomalousYX = py/mass*elec_charge*d0*ny*D/permitivity;

    return px * py / mass - odd_vis*Vxx + anomalousYX;
}

float BerryFluid::YMomentumFluxY(StateVec2D U) {
    float sound=U.S();
    float den=U.n();
    float py=U.py();
    float mass=DensityToMass(den);
    float Vyx=U.dyvx();
    float nx=U.dxn();
    float temp = U.tmp();
    //float D = AnomalousStressTensor(temp, chem_pot, gap);
    float D = 0.1;
    float anomalousYY = -py/mass*elec_charge*d0*nx*D/permitivity;

    return py * py / mass + vel_fer * vel_fer * mass / 3.0f + 0.5f * sound * sound * den * den - odd_vis*Vyx + anomalousYY;
}

BerryFluid::~BerryFluid(){
delete[] Umain;
delete[] Umid;
delete[] Den;
delete[] VelX;
delete[] VelY;
delete[] CurX;
delete[] CurY;
delete[] Tmp;
}


