
/*!@file
 * @brief Header file for 1D Dirac-Fermi fluid
 */


#ifndef DIRACGRAPHENE1DLIB_H
#define DIRACGRAPHENE1DLIB_H

#include "includes/Fluid1DLib.h"


/*!
 * @brief Graphene electronic fluid class in one dimension.
 *
 * The DiracGraphene1D class describes ...
 * It overrides class Fluid1D necessary methods in order to describe the semi-classical electronic fluid.
 * */
class DiracGraphene1D : public Fluid1D{
protected:
    /*!
	 * @brief Converts the number density to efective mass
	 *
	 * Since the mass of the fluid element is not a constant in the graphene electronic fluid, one needs to perform the transformation
	   @f[ m^\star = n^{3/2} ]
	 * */
    float DensityToMass(float density) override;


    //StateVec1D *ptr_StateVecHole = nullptr;


    float vel_therm = 10.0f; //new constant - pressure term
    float A = 0.1f; //new constant - source function, equilibrium relaxation
    float B = 1.0f; //new constant - source function, electron-hole creation

public:

    StateVec1D * HoleUmain;
    StateVec1D * HoleUmid;

    Group* GrpHDen ;    ///< group for ALL Hole Density snapshots
    Group* GrpHVel ;    ///< group for ALL Hole Velocity snapshots
    DataSpace* DataspaceHDen;    ///< dataspace for EACH Hole Density snapshots
    DataSpace* DataspaceHVel;    ///< dataspace for EACH Hole Velocity snapshots

    float *HDen;    // hole number density
    float *HVel;    // hole fluid velocity

    explicit DiracGraphene1D(SetUpParametersCNP &input_parameters);
        ~DiracGraphene1D();

        void SetSound() override;
        void CflCondition() override;
        void InitialCondRand() override;
        //void SetSimulationTime() override;


    void RichtmyerStep1() override;
    void RichtmyerStep2() override;


        float EleDensitySource(StateVec1D Uelec , StateVec1D Uholes);   ///< density equation (continuity equation) source term
        float EleVelocitySource(StateVec1D Uelec); ///< velocity X component equation (momentum equation) source term

        float EleDensityFlux(StateVec1D Uelec); ///< density equation (continuity equation) conserved flux
        float EleVelocityFlux(StateVec1D Uelec); ///< velocity X component equation (momentum equation) conserved flux

        float HolDensitySource(StateVec1D Uelec , StateVec1D Uholes);   ///< density equation (continuity equation) source term
        float HolVelocitySource(StateVec1D Uholes); ///< velocity equation (momentum equation) source term

        float HolDensityFlux(StateVec1D Uholes); ///< density equation (continuity equation) conserved flux
        float HolVelocityFlux(StateVec1D Uholes); ///< velocity equation (momentum equation) conserved flux

        void WriteFluidFile(float t) override; // writes the line of time t on the simplified .dat file output

        //void ForwardTimeOperator();

    void CopyFields() override;
    void SaveSnapShot() override;
    void CreateHdf5File() override;

};

#endif //DIRACGRAPHENE1DLIB_H