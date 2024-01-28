//
// Created by pedro on 07-01-2024.
//

#ifndef TRIPLEGRAPHENE2DLIBNOT_H
#define TRIPLEGRAPHENE2DLIBNOT_H

#include "includes/DiracGraphene2DLib.h"

class TripleGraphene2D : public DiracGraphene2D {
protected:

    StateVec2D *ptr_StateVecBoson = nullptr;

    void ChooseGridPointers(const string &grid) override;

public:
    StateVec2D * BUmain;
    StateVec2D * BUmid;

    Group* GrpBDen ;     ///< group for ALL Boson Density snapshots
    Group* GrpBVelX ;    ///< group for ALL Boson Velocity X snapshots
    Group* GrpBVelY ;    ///< group for ALL Boson Velocity X snapshots
    DataSpace* DataspaceBDen;    ///< dataspace for EACH Boson Density snapshots
    DataSpace* DataspaceBVelX;   ///< dataspace for EACH Boson Velocity X snapshots
    DataSpace* DataspaceBVelY;   ///< dataspace for EACH Boson Velocity Y snapshots

    float *BDen;       // boson number density
    float *BVelX;      // boson fluid velocity x component
    float *BVelY;      // boson fluid velocity y component


    explicit TripleGraphene2D(SetUpParametersCNP &input_parameters);
    ~TripleGraphene2D();

    void SetSound() override;
    float FermionDensityToMass(float density); ///<calculo da massa hidrodinamica para fermioes
    float BosonDensityToMass(float density); ///< calculo da massa hidrodinamica para bosoes
    void InitialCondRand() override; ///< condicao inicial de densidades aleatorias
    void InitialCondPointDen() override; ///< permite establecer condicoes inicias com os indices da grelha de simulacao
    void InitialCondTest() override; ///< permite estabelecer condicoes inicias em coordenadas cartesianas para a densidade

    void RichtmyerStep1() override;
    void RichtmyerStep2() override;

    float EleDensitySource(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson);   ///< density equation (continuity equation) source term
    float EleXMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson); ///< velocity X component equation (momentum equation) source term
    float EleYMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity y component equation (momentum equation) source term

    float EleDensityFluxX(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,StateVec2D Uboson); ///< density equation (continuity equation) conserved flux X component
    float EleDensityFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson) ; ///< density equation (continuity equation) conserved1 flux Y component

    float EleXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux X component
    float EleXMomentumFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux Y component

    float EleYMomentumFluxX(StateVec2D Uelec , __attribute__((unused))  StateVec2D Uholes,  StateVec2D Uboson); ///< velocity Y component equation (momentum equation) conserved flux X component
    float EleYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) ; ///< velocity Y component equation (momentum equation) conserved flux Y component

    float HolDensitySource(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson);   ///< density equation (continuity equation) source term
    float HolXMomentumSource(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity X component equation (momentum equation) source term
    float HolYMomentumSource(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity y component equation (momentum equation) source term

    float HolDensityFluxX(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson); ///< density equation (continuity equation) conserved flux X component
    float HolDensityFluxY(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes,  StateVec2D Uboson); ///< density equation (continuity equation) conserved1 flux Y component

    float HolXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes,StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux X component
    float HolXMomentumFluxY(__attribute__((unused))  StateVec2D Uelec , StateVec2D Uholes,StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux Y component

    float HolYMomentumFluxX(__attribute__((unused)) StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson); ///< velocity Y component equation (momentum equation) conserved flux X component
    float HolYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson); ///< velocity Y component equation (momentum equation) conserved flux Y component

    float BosonDensitySource(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson);   ///< density equation (continuity equation) source term
    float BosonXMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson); ///< velocity X component equation (momentum equation) source term
    float BosonYMomentumSource( __attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity y component equation (momentum equation) source term

    float BosonDensityFluxX(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes, StateVec2D Uboson); ///< density equation (continuity equation) conserved flux X component
    float BosonDensityFluxY(__attribute__((unused)) StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson) ; ///< density equation (continuity equation) conserved1 flux Y component

    float BosonXMomentumFluxX(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux X component
    float BosonXMomentumFluxY(StateVec2D Uelec , __attribute__((unused)) StateVec2D Uholes,  StateVec2D Uboson); ///< velocity X component equation (momentum equation) conserved flux Y component

    float BosonYMomentumFluxX(StateVec2D Uelec , __attribute__((unused))  StateVec2D Uholes,  StateVec2D Uboson); ///< velocity Y component equation (momentum equation) conserved flux X component
    float BosonYMomentumFluxY(StateVec2D Uelec , StateVec2D Uholes, StateVec2D Uboson) ; ///< velocity Y component equation (momentum equation) conserved flux Y component

    void WriteFluidFile(float t) override; ///<  writes the line of time t on the simplified .dat file output

    void ForwardTimeOperator() override; ///< Time evolution for the FTCS method employed for the parabolic operators.
    void ForwardTimeOperator(char field) override; ///< Time evolution for the FTCS method employed for the parabolic operators.

    void CopyFields() override;
    void SaveSnapShot() override;
    void CreateHdf5File() override;


};


#endif //TRIPLEGRAPHENE2DLIBNOT_H
