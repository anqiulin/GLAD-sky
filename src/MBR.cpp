/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/
#include "MBR.h"


MBR::MBR()
{
}

MBR::MBR(int _minx, int _miny, int _maxx, int _maxy)
{
  minx = _minx; miny = _miny; maxx = _maxx; maxy = _maxy;
}

MBR::MBR(Point& minPoint, Point& maxPoint)
{
  minx = minPoint.coord[0]; miny = minPoint.coord[1]; maxx = maxPoint.coord[0]; maxy = maxPoint.coord[1];
}
void MBR::unionWith(Point* p)
{
  minx = MIN(minx, p->coord[0]);
  maxx = MAX(maxx, p->coord[0]);
  miny = MIN(miny, p->coord[1]);
  maxy = MAX(maxy, p->coord[1]);

}
int MBR::area()
{
  return (maxx - minx) * (maxy - miny);
}
int MBR::width()
{
  return maxx - minx;
}
int MBR::height()
{
  return maxy - miny;
}
double MBR::diagonal()
{
  return sqrt((maxy - miny)*(maxy - miny) + (maxx - minx)*(maxx - minx));
}

MBR::~MBR()
{
}
