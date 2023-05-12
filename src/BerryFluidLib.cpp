/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281                                                                    *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


#include "includes/BerryFluidLib.h"

float theta = .0;
float A = .01;

BerryFluid::BerryFluid(SetUpParameters &input_parameters) : Fluid2D(input_parameters) {
    //vel_fer = input_parameters.FermiVelocity ;//fermi_velocity;
    //therm_diff = input_parameters.ThermalDiffusivity; //thermal diffusivity
    //theta = input_parameters.ElectricFieldAngle; //angle described by the electric field in relation to xx
    //delta = input_parameters.BandGap; //band gap
    //miu = input_parameters.ChemicalPotential
    char buffer [100];
    sprintf (buffer, "S=%.2fvF=%.2fvis=%.3fodd=%.3fl=%.3fwc=%.2ftherm=%.2f", vel_snd, vel_fer, kin_vis,odd_vis, col_freq,cyc_freq,therm_diff);
    file_infix = buffer;
}


float BerryFluid::XMomentumFluxX(StateVec2D U) {
    float den=U.n(); //aqui não seria para usar a relação n=1/(2*M_PI*hbar^2*v_F^2)*(miu^2-delta^2)? e nesse caso o miu é a energia de fermi?
    float px=U.px();

    return U.px()*U.px()/U.n() + U.n()+A*px/pow(den, 3/2)*sin(theta);
}

float BerryFluid::YMomentumFluxX(StateVec2D U) {
    float den=U.n();
    float py=U.py();

    return U.py()*U.px()/U.n() + A*py/pow(den, 3/2)*sin(theta);
}

float BerryFluid::XMomentumFluxY(StateVec2D U) {
    float den=U.n();
    float px=U.px();

    return U.px()*U.py()/U.n() - A*px/pow(den, 3/2)*cos(theta);
}

float BerryFluid::YMomentumFluxY(StateVec2D U) {
    float den=U.n();
    float py=U.py();

    return U.py()*U.py()/U.n() + U.n() - A*py/pow(den, 3/2)*cos(theta);
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