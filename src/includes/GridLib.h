
#ifndef GRIDLIB_H
#define GRIDLIB_H

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
	int NW;
	int NE;
	int SW;
	int SE;
	int N;
	int S;
	int E;
	int W;
	int C;
	void FirstNeighbours();
	void SecondNeighbours();
};


#endif //GRIDLIB_H
