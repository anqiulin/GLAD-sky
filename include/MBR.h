/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/

#ifndef MBR_H
#define MBR_H
#include"Point.h"

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ( (x) > (y) ? (x) : (y) )

class MBR
{
public:
	int minx = 999999999;
	int miny = 999999999;
	int maxx = -999999999;
	int maxy = -999999999;
	MBR();
	MBR(int _minx, int _miny, int _maxx, int _maxy);
	MBR(Point & minPoint, Point & maxPoint);
	void unionWith(Point* _point);
	int area();
	int width();
	int height();
	double diagonal();
	~MBR();
};
#endif
