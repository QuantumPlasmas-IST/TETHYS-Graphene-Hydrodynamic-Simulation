/************************************************************************************************\
* 2020 Pedro Cosme , João Santos and Ivan Figueiredo                                             *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#ifndef GRID2DLIB_H
#define GRID2DLIB_H

#include <cstdlib>

class GridPoint{
private:
	int size_x;
	int size_y;
	int stride;

public:
	GridPoint(int pos,int xpoints, int ypoints,bool mid);
	~GridPoint() = default;
	bool IsMidGrid = false;
	int NW=0; //primeiros vizinhos portanto pertencentes a outro tipo de grelha
	int NE=0;
	int SW=0;
	int SE=0;
	int NW2=0; //terceiros vizinhos mesmo tipo de grelha
	int NE2=0;
	int SW2=0;
	int SE2=0;
	int N=0; //segundos vizinhos mesmo tipo de grelha
	int S=0;
	int E=0;
	int W=0;
	int C=0; //ponto de referencia da grelha
	void FirstNeighbours();
	void SecondNeighbours();
	void ThirdNeighbours();
};


#endif //GRID2DLIB_H
