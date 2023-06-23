//
// Created by pedro on 05-05-2023.
//

#ifndef TRIPLEGRAPHENE1DLIB_H
#define TRIPLEGRAPHENE1DLIB_H


#include "includes/DiracGraphene1DLib.h"

/*!
 * @brief Graphene electronic fluid class in one dimension.
 *
 * The TripleGraphene1D class describes ...
 * It overrides class FLuid1D and DiracGraphene1D necessary methods in order to describe the semi-classical electronic fluid.
 * */


class TripleGraphene1D : public DiracGraphene1D{
protected:

    float vel_therm = 10.0f; //new constant - pressure term
    float A = 0.1f; //new constant - source function, equilibrium relaxation
    float B = 1.0f; //new constant - source function, electron-hole creation

public:

    StateVec1D * BUmain;
    StateVec1D * BUmid;

    Group* GrpBDen ;    ///< group for ALL Hole Density snapshots
    Group* GrpBVel ;    ///< group for ALL Hole Velocity snapshots
    DataSpace* DataspaceBDen;    ///< dataspace for EACH Hole Density snapshots
    DataSpace* DataspaceBVel;    ///< dataspace for EACH Hole Velocity snapshots

    float *BDen;    // hole number density
    float *BVel;    // hole fluid velocity

    explicit TripleGraphene1D(SetUpParametersCNP &input_parameters);
        ~TripleGraphene1D();

    void SetSound() override;
    //void CflCondition() override;
    void InitialCondRand() override;
    //void SetSimulationTime() override;


    //void ComputeElectricPotencial(const string &grid);
    //float PotencialKernel(float x);


    void RichtmyerStep1() override;
    void RichtmyerStep2() override;


    float EleDensitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes);   ///< density equation (continuity equation) source term
    float EleVelocitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< velocity X component equation (momentum equation) source term

    float EleDensityFlux(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< density equation (continuity equation) conserved flux
    float EleVelocityFlux(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< velocity X component equation (momentum equation) conserved flux

    float HolDensitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes);  ///< density equation (continuity equation) source term
    float HolVelocitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< velocity equation (momentum equation) source term

    float HolDensityFlux(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< density equation (continuity equation) conserved flux
    float HolVelocityFlux(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< velocity equation (momentum equation) conserved flux

    float BDensitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes);  ///< density equation (continuity equation) source term
    float BVelocitySource(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes);///< velocity X component equation (momentum equation) source term

    float BDensityFlux(StateVec1D Uboson , StateVec1D Uelec, StateVec1D Uholes); ///< density equation (continuity equation) conserved flux
    float BVelocityFlux(StateVec1D Uboson, StateVec1D Uelec, StateVec1D Uholes); ///< velocity equation (momentum equation) conserved flux

    void WriteFluidFile(float t) override; // writes the line of time t on the simplified .dat file output

    //void ForwardTimeOperator();

    void CopyFields() override;
    void SaveSnapShot() override;
    void CreateHdf5File() override;

    void InitialCondTest() override;


};



#endif //TRIPLEGRAPHENE1DLIB_H
