/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/
#ifndef GRID_H
#define GRID_H
#include<map>
#include<vector>
#include<set>
#include<string.h>
#include<math.h>
#include"Point.h"
#include"MBR.h"
#include"Graph.h"
#include<algorithm>

using namespace std;
class Grid
{
public:
	Grid(int _colSize, int _rowSize);
	Grid(int _colSize, int _rowSize, MBR* _extend);
	Grid(MBR* _extend);
	Grid();
	Grid(Grid* grid);
	bool isBuilt();
	int numCells();
	~Grid();
	vector<vector<int>> nodeMap;
	struct gridNode{
		bitset<KIND> gridNodeBits;
		vector<int> POI;
	};
	vector<gridNode> gridNodeMap;
	int colSize;	// # of columns
	int rowSize;	// # of rows
	int colLength;
	int rowLength;
	MBR *extend;
	int* cellContain; // a tag to check whether a cell contains objects
	void getcell_Straight(Point *s, Point *t, bitset<KIND> query, vector<int>* result);
	void getCell_MBR(Point *s, Point *t, bitset<KIND> cellbit, vector<int>* result);
	void getPoiInCells(vector <int>* cell, vector<int>* node);
	int getCell(Point *p); // compute the cell where the given point belongs to
	set<int> getNeighbors(int cid, set<int> result);
	set<int> getNeighborsMultiRound(int cid, set<int> rst, int round);
	void getNodeInCells(vector<int>* _cells, vector<int>* result, int & objectNum);
	void getKObjectNode(int _k, int _cid, vector<int>* result);
	void getObjectInCells(vector<int>* _cells, vector<int>* result);
	void getKObjects(int _k, int _cid, vector<int>* result);
	set<int> getRoundCellsSet(int _round, int _cid);
	void getRoundCells(int _round, int _cid, vector<int>* result);
	void calculateGridSize(int decimeter);
	void build(Point** _pts, Graph* G, int decimeter);

	int rightOf(int cid); // return the right cell of given cell
	int upperOf(int cid); // return the upper cell of given cell
	int belowOf(int cid); // return the below cell of given cell
	int leftOf(int cid); // return the left cell of given cell
};
#endif
