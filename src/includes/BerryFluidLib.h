#ifndef BERRYFLUIDLIB_H
#define BERRYFLUIDLIB_H

#include "includes/Fluid2DLib.h"


class BerryFluid : public Fluid2D{
protected:
    float DensityToMass(float density)override;

public :
    explicit BerryFluid(SetUpParameters &input_parameters);
        ~BerryFluid();

	void CflCondition() override;

    float XMomentumSource(StateVec2D U)override; ///< velocity X component equation (momentum equation) source term
    float YMomentumSource(StateVec2D U)override; ///< velocity y component equation (momentum equation) source term

    float XMomentumFluxX(StateVec2D U )override; ///< velocity X component equation (momentum equation) conserved flux X component
    float XMomentumFluxY(StateVec2D U )override; ///< velocity X component equation (momentum equation) conserved flux Y component

    float YMomentumFluxX(StateVec2D U)override; ///< velocity Y component equation (momentum equation) conserved flux X component
    float YMomentumFluxY(StateVec2D U )override; ///< velocity Y component equation (momentum equation) conserved flux Y component

};


#endif //BERRYFLUIDLIB_H
