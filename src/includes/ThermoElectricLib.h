/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

//
// Created by pcosme on 22/10/2021.
//

/*!@file
 * @brief Header file for dealing with thermo electric effects
 */

#ifndef THERMOELECTRICLIB_H
#define THERMOELECTRICLIB_H

#include "includes/TethysBaseLib.h"

class ThermoElectric {
	public:
	static float ThermoPower(float density,float temperature);
};


#endif //THERMOELECTRICLIB_H
