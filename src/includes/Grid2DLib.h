/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/


/*!@file
 * @brief Header file for grid point class
 */

#ifndef GRID2DLIB_H
#define GRID2DLIB_H

#include <cstdlib>

/*!
 * @brief Ancillary class for dealing with the mesh points.
 *
 * This class is responsible for finding the position of the grid neighbours, either at the main or mid grid, for the numerical methods employed.
 * */
class GridPoint{
private:
	int size_x;
	int size_y;
	int stride;

public:
	GridPoint(int pos,int xpoints, int ypoints,bool mid);
	~GridPoint() = default;
	bool IsMidGrid = false;
	int NW=0; //second neighbours -> different type of grid
	int NE=0;
	int SW=0;
	int SE=0;
	int NW2=0; //third neighbours -> same type of grid
	int NE2=0;
	int SW2=0;
	int SE2=0;
	int N=0; //second neighbours -> same type of grid
	int S=0;
	int E=0;
	int W=0;
	int C=0; //central point
	void FirstNeighbours(); ///< Finds the first neighbours of a mesh point, i.e. sets the NW NE SW and SE points
	void SecondNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the N S E and W points
	void ThirdNeighbours(); ///< Finds the second neighbours of a mesh point, i.e. sets the NW2 NE2 SW2 and SE2 points
};


#endif //GRID2DLIB_H
