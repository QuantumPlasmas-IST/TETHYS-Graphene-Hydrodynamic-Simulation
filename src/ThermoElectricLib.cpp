/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

//
// Created by pcosme on 22/10/2021.
//

#include "ThermoElectricLib.h"

float ThermoElectric::ThermoPower(float density, float temperature) {
	return -PHYS_SEEBECK_0*temperature/sqrt(density);
}
