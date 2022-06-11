/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for initialization class
 */

#ifndef SETUPPARAMETERSLIB_H
#define SETUPPARAMETERSLIB_H
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <string>
#include <random>
#include <regex>
#include <exception>
#include <H5Cpp.h>
#include <omp.h>

/*!
 * @brief Initialization class for the fluid classes.
 *
 * This class allows the initialization of the appropriate physical parameters either from inline arguments, prompt or reading them from an existing HDF5 file
 * */
class SetUpParameters {
	public:
		SetUpParameters();
		SetUpParameters(int argc, char ** argv);
		SetUpParameters(float sound, float fermi, float coll, float visco, float odd_visco, float cyclo,
		                float thermal_diff, int mode, float aspect);
		~SetUpParameters() = default;
		int SaveMode;
		int SizeX;
		int SizeY;
		float Length=1.0f;
		float Width=1.0f;
		float AspectRatio=1.0f;
		float SoundVelocity;
		float FermiVelocity;
		float CollisionFrequency;
		float ShearViscosity;
		float OddViscosity;
		float CyclotronFrequency;
		float ThermalDiffusivity;
	    float SimulationTime=0.0f;
		void ParametersChecking() const; ///< Runs a checking on the physical feasibility of the parameters
		/*!
		 * @brief Sets up the 2D grid dimensions
		 *
		 * The rectangular simulation domain is characterized by the aspect ratio @f$ AR=L/W @f$ and the number of points of the smallest dimension @f$ N @f$.
		 * @f[ AR\geq1 \Rightarrow N_x=AR(N_y-1)+1 \quad N_y=N @f]
		 * @f[ AR<1 \Rightarrow N_x=N \quad N_y=\frac{N_x-1}{AR}+1 @f]
		 * */
		void DefineGeometry();
		void ParametersFromHdf5File(const std::string& hdf5name); ///< Imports the parameters from a saved HDF5 file
        void PrintParameters() const; ///< Prints the read parameters to standard output
		void PromptParameters() ; ///< Asks the user for the simulation parameters
		/*!
		 * @brief Imports the parameters from a saved .ini file
		 *
		 * As standard with .ini files sections can be indicated by square brackets as in `[section]` and comment lines start with a semicolon `;`
		 * The readble keywords are
		 *
		 *
		 * | Keywords | Parameter |   |
		 * | ----: | :----: | :---: |
		 * |  sound   | Sound velocity     | @f$ S @f$    |
		 * |  fermi   | Fermi velocity     | @f$ v_0 @f$    |
		 * |  shear   | Shear viscosity     | @f$ \nu_s @f$    |
		 * |  odd   |  Odd/Hall viscosity    | @f$ \nu_o @f$    |
		 * |  col   |  Collision frequency    | @f$ 1/\tau @f$    |
		 * |  cycl   |  Cyclotron frequency    | @f$ \omega_c @f$    |
		 * |  therm   |  Thermal diffusivity    | @f$ \alpha @f$    |
		 * |  aspect   |  Aspect ratio    | @f$ AR @f$    |
		 * |  save   |  Save mode    | -    |
		 * */
        void ReadIniFile(char * file_name);
};



#endif //SETUPPARAMETERSLIB_H
