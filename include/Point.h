/*
	Author: Qiulin An
	Date: 6/5/2023
*/

#ifndef POINT_H
#define POINT_H
#include<math.h>
#include<cmath>

class Point
{
public:
	int id;
	int coord[2];
	int nodeId = -1;
	Point();
	Point(int _x, int _y);
	Point(int _id, int _x, int _y);
	double Dist(Point _point);
	double deg2rad(double deg);
	double rad2deg(double rad);
	double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d);
	~Point();
};
#endif
