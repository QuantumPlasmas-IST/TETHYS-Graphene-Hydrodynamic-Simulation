/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for grid point class
 */

#ifndef GRIDLIB_H
#define GRIDLIB_H

#include <cstdlib>

/*!
 * @brief Ancillary class for dealing with the mesh points.
 *
 * This class is responsible for finding the position of the grid neighbours, either at the main or mid grid, for the numerical methods employed.
 * */
	class GridPoint2D {
	private:
		int size_x;
		int size_y;
		int stride;

	public:
		GridPoint2D(int pos, int xpoints, int ypoints, bool mid);

		~GridPoint2D() = default;

		bool IsMidGrid = false;
		int NW = 0; //first neighbours -> different type of grid
		int NE = 0;
		int SW = 0;
		int SE = 0;
		int NW2 = 0; //third neighbours -> same type of grid
		int NE2 = 0;
		int SW2 = 0;
		int SE2 = 0;
		int N = 0; //second neighbours -> same type of grid
		int S = 0;
		int E = 0;
		int W = 0;
		int C = 0; //central point
		void FirstNeighbours(); ///< Finds the first neighbours of a mesh point, i.e. sets the NW NE SW and SE points
		void SecondNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the N S E and W points
		void
		ThirdNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the NW2 NE2 SW2 and SE2 points
	};



/*!
 * @brief Ancillary class for dealing with the mesh points.
 *
 * This class is responsible for finding the position of the grid neighbours, either at the main or mid grid, for the numerical methods employed.
 * */
	class GridPoint1D {
	private:
		int size_x;

	public:

		GridPoint1D(int pos, int xpoints, bool mid);
		~GridPoint1D() = default;

		bool IsMidGrid = false;
		int C = 0; //central point
		int W = 0; //first neighbours -> different type of grid
		int E = 0;
		int W2 = 0; //secon neighbours -> same type of grid
		int E2 = 0;
		int W3 = 0; //third neighbours -> same type of grid
		int E3 = 0;

		void FirstNeighbours(); ///< Finds the first neighbours of a mesh point, i.e. sets the W and E points
		void SecondNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the WW and EE points
		void ThirdNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the WWW and EEE points
	};



#endif //GRIDLIB_H
