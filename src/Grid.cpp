/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/
#include "Grid.h"
#include "Graph.h"
#include <sys/time.h>


static const int MINPOINT = 4;


Grid::Grid(int _colSize, int _rowSize)
{
	colSize = _colSize;
	rowSize = _rowSize;
}


Grid::Grid(Grid* grid){
	colSize = grid->colSize;
	rowSize = grid->rowSize;
	//cellSize = grid->cellSize;
	extend = grid->extend;
	nodeMap = grid->nodeMap;
	cellContain = new int[colSize * rowSize];
	memcpy(cellContain, grid->cellContain, sizeof(grid->cellContain));
}

Grid::Grid()
{
}

Grid::Grid(MBR* _extend)
{
	extend = _extend;
}

bool Grid::isBuilt() {
	return (nodeMap.size()>0);
}

int Grid::numCells() {
	return rowSize*colSize;
};

int Grid::getCell(Point *p)
{
	int col = (p->coord[0] - extend->minx) / colLength;
	int row = (p->coord[1] - extend->miny) / rowLength;
	return  row*colSize + col;
}
void Grid::getcell_Straight(Point *s, Point *t, bitset<KIND> query, vector<int>* result)
{
	int dx = t->coord[0] - s->coord[0];
	int dy = t->coord[1] - s->coord[1];
	vector<int> tmp;
	
	int increx = dx/abs(dx);
	int increy = dy/abs(dy);

	double k = (dy*1.0)/(dx*1.0);
	double b = -k * (s->coord[0] * 1.0) + (s->coord[1] *1.0);
	// cout << dx << endl;
	bitset<KIND>cellbit;
	for(int i = 0; i < abs(dx); i += colLength/3)
	{
		int x = s->coord[0] + increx*i;
		int y = (int)(k * x + b); 

		int col = (x - extend->minx) / colLength;
		int row = (y - extend->miny) / rowLength;

		int cell = row*colSize + col;
		tmp.push_back(cell);

		if(cellContain[cell] > 0)
		{
			if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
			{
				cellbit |= gridNodeMap[cell].gridNodeBits & query;
				result->push_back(cell);
			}	
		}
	}
	for(int i = 0; i <= abs(dy); i += rowLength/3)
	{
		int y = s->coord[1] + increy*i;
		int x = (int)((y - s->coord[1])/k) + s->coord[0];


		int col = (x - extend->minx) / colLength;
		int row = (y - extend->miny) / rowLength;

		int cell = row*colSize + col;
		tmp.push_back(cell);

		if(cellContain[cell] > 0)
		{
			if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
			{
				cellbit |= gridNodeMap[cell].gridNodeBits & query;
				result->push_back(cell);
			}	
		}
	}
	sort(tmp.begin(),tmp.end());
	tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
	vector<int>temp;
	// for(auto c : tmp)
	// {
	// 	cout << c << " ";
	// }
	// cout << endl;
	// cout << cellbit.count() << " " << query.count() << endl;
	while(cellbit.count() != query.count())
	{
		// cout << cellbit.count() << " " << query.count() << endl;
		temp.clear();
		for(int i = 0; i < tmp.size(); i++)
		{
			int cell = belowOf(tmp[i]);
			if(cell != -1)
			{
				temp.push_back(cell);
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}	
				}
			}
			cell = upperOf(tmp[i]);
			if(cell != -1)
			{
				temp.push_back(cell);
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}	
				}
			}
			cell = leftOf(tmp[i]);
			if(cell != -1)
			{
				temp.push_back(cell);
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}	
				}
			}
			cell = rightOf(tmp[i]);
			if(cell != -1)
			{
				temp.push_back(cell);
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}	
				}
			}
		}
		sort(temp.begin(),temp.end());
		temp.erase(unique(temp.begin(), temp.end()), temp.end());
		// for(auto c : temp)
		// {
		// 	cout << c << " ";
		// }
		// cout << endl;
		tmp.clear();
		tmp.insert(tmp.begin(),temp.begin(),temp.end());
	}
	sort(result->begin(),result->end());
	result->erase(unique(result->begin(), result->end()), result->end());
}
void Grid::getCell_MBR(Point *s, Point *t, bitset<KIND> query, vector<int>* result)
{
	bitset<KIND> cellbit;
	cellbit.reset();
	result->clear();

	int s_col = (s->coord[0] - extend->minx) / colLength;
	int s_row = (s->coord[1] - extend->miny) / rowLength;

	int t_col = (t->coord[0] - extend->minx) / colLength;
	int t_row = (t->coord[1] - extend->miny) / rowLength;

	int low_col = min(s_col, t_col);
	int low_row = min(s_row, t_row);

	int up_col = max(s_col, t_col);
	int up_row = max(s_row, t_row);

	for(int i = low_row; i <= up_row; i++)
	{
		for(int j = low_col; j <= up_col; j++)
		{
			int cell = i*colSize + j;
			if(cellContain[cell] > 0)
			{
				if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
				{
					cellbit |= gridNodeMap[cell].gridNodeBits & query;
					result->push_back(cell);
				}
				
			}
		}
	}
	int flag = 1;
	while(flag || query.count() != cellbit.count())
	{
		flag = 0;
		int round = 1;
		low_col -= round;
		low_row -= round;
		up_col += round;
		up_row += round;

		if(low_col >= 0)
		{
			for(int i = low_row; i <= up_row; i++)
			{
				if(i < 0)
					continue;
				if(i >= rowSize)
					break;
				
				int cell = i*colSize + low_col;
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}
				}
			}
		}
		if(up_col < colSize)
		{
			for(int i = low_row; i <= up_row; i++)
			{
				if(i < 0)
					continue;
				if(i >= rowSize)
					break;
				
				int cell = i*colSize + up_col;
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}
				}
			}
		}
		if(low_row >=0)
		{
			for(int i = low_col+1; i <= up_col-1; i++)
			{
				if(i < 0)
					continue;
				if(i >= colSize)
					break;
				
				int cell = low_row *colSize + i;
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}
				}
				
			}
		}
		if(up_row < rowSize)
		{
			for(int i = low_col+1; i <= up_col-1; i++)
			{
				if(i < 0)
					continue;
				if(i >= colSize)
					break;
				
				int cell = up_row *colSize + i;
				if(cellContain[cell] > 0)
				{
					if((gridNodeMap[cell].gridNodeBits & query).count() > 0)
					{
						cellbit |= gridNodeMap[cell].gridNodeBits & query;
						result->push_back(cell);
					}
				}
				
			}
		}
	}
}

void Grid::calculateGridSize(int decimeter)
{
	int minx, miny, maxx, maxy;
	minx = extend->minx;
	miny = extend->miny;
	maxx = extend->maxx;
	maxy = extend->maxy;
	Point leftDown;
	leftDown.coord[0] = minx;
	leftDown.coord[1] = miny;
	Point leftUp;
	leftUp.coord[0] = minx;
	leftUp.coord[1] = maxy;
	Point rightDown;
	rightDown.coord[0] = maxx;
	rightDown.coord[1] = miny;
	double xLength = leftDown.Dist(rightDown);
	double yLength = leftDown.Dist(leftUp);
	colLength = (maxx - minx) * ((double)decimeter/xLength);
	rowLength = (maxy - miny) * ((double)decimeter/yLength);
	cout << "minx = " << minx << "\n";
	cout << "maxx = " << maxx << "\n";
	cout << "miny = " << miny << "\n";
	cout << "maxy = " << maxy << "\n";
	cout << "xLength = " << xLength << "\n";
	cout << "yLength = " << yLength << "\n";
	cout << "colLength = " << colLength << "\n";
	cout << "rowLength = " << rowLength << "\n";
	cout << "decimeter = " << decimeter << "\n";
}

void Grid::build(Point** _pts, Graph* G, int decimeter)
{
	if (G->POIObjects.size() < MINPOINT)
		return;
	calculateGridSize(decimeter);
	colSize = ceil((double)extend->width() / colLength);
	cout<<"colSize = " << colSize <<endl;
	rowSize = ceil((double)extend->height() / rowLength);
	cout<< "rowSize = " << rowSize << endl;
	int numNode = colSize * rowSize;
	cout<< "numNode = " << numNode << endl;
	cellContain = new int[numNode];
	
	for(int i = 0; i < numNode; i++)
	{
		gridNode item;
		gridNodeMap.push_back(item);
		cellContain[i] = 0;
	}
	
	for (int pid = 0; pid < G->POIObjects.size(); pid++)
	{
		int gridid = getCell(_pts[pid]);
		int node = _pts[pid]->nodeId;

		gridNodeMap[gridid].gridNodeBits |= G->NodesBit[node];
		gridNodeMap[gridid].POI.push_back(node);
		cellContain[gridid]++;

		if(gridid >= numNode)
			printf("Error: Grid::build() occurs an error!\n");
	}
}

void Grid::getPoiInCells(vector <int>* cell, vector<int>* node)
{
	for(int i = 0; i < cell->size(); i++)
	{
		node->insert(node->begin(), gridNodeMap[cell->at(i)].POI.begin(), gridNodeMap[cell->at(i)].POI.end());
	}
}

// void Grid::getKObjectNode (int _k, int _cid, vector<int>* result)
// {
// 	if(_cid >= colSize * rowSize)
// 		return;
// 	int round = 0;
// 	int objectNum = 0;
// 	if(cellContain[_cid] > 0){
// 		for(int i = 0; i < cellContain[_cid]; i++){
// 			result->push_back(gridNodeMap[_cid][i].nodeId);
// 			objectNum = objectNum + gridNodeMap[_cid][i].objects.size();
// 		}
// 	}
// 	while(objectNum < _k)
// 	{
// 		round++;
// 		vector<int> cells;
// 		getRoundCells(round, _cid, &cells);

// 		vector<int> nodes;
// 		getNodeInCells(&cells, &nodes, objectNum);
// 		if(nodes.size() > 0) {
// 			result->insert(result->end(), nodes.begin(), nodes.end());
// 		}
// 		cells.clear();
// 		nodes.clear();
// 	}
// 	result->push_back(round);
// }

void Grid::getKObjects (int _k, int _cid, vector<int>* result)
{
	if(_cid >= colSize * rowSize)
		return;
	int round = 0;
	if(cellContain[_cid] > 0){
		result->insert(result->end(), nodeMap[_cid].begin(), nodeMap[_cid].end());
	}
	while(result->size() < _k)
	{
		round++;
		vector<int> cells;
		getRoundCells(round, _cid, &cells);
		vector<int> objects;
		getObjectInCells(&cells, &objects);
		if(objects.size() > 0) {
			result->insert(result->end(), objects.begin(), objects.end());
		}
		cells.clear();
		objects.clear();
	}
	result->push_back(round);
}


void Grid::getRoundCells(int _round, int _cid, vector<int>* result)
{
	if(_cid >= colSize * rowSize)
		return;
	int x = _cid % colSize; // column
	int y = _cid / colSize; // row
	int leftx = x - _round;
	int lowy = y - _round;
	int rightx = x + _round;
	int upy = y + _round;
	if(leftx >= 0)
	{
		int cell = y * colSize + leftx;
		if(cellContain[cell] > 0)
			result->push_back(cell);
	}
	if(rightx < colSize)
	{
		int cell = y * colSize + rightx;
		if(cellContain[cell] > 0)
			result->push_back(cell);
	}
	if(lowy >= 0)
	{
		int cell = lowy * colSize + x;
		if(cellContain[cell] > 0)
			result->push_back(cell);
	}
	if(upy < rowSize)
	{
		int cell = upy * colSize + x;
		if(cellContain[cell] > 0)
			result->push_back(cell);
	}
	for(int i = 0; i < _round; i++)
	{
		if(leftx >= 0)
		{
			if((y - (i + 1)) >= 0){
				int cell = (y - (i + 1)) * colSize + leftx;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
			if((y + (i + 1)) < rowSize){
				int cell = (y + (i + 1)) * colSize + leftx;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
		}
		if(rightx < colSize)
		{
			if((y - (i + 1)) >= 0){
				int cell = (y - (i + 1)) * colSize + rightx;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
			if((y + (i + 1)) < rowSize){
				int cell = (y + (i + 1)) * colSize + rightx;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
		}
	}
	for(int i = 1; i < _round; i++)
	{
		if(lowy >= 0)
		{
			if((x - i) >= 0){
				int cell = lowy * colSize + x - i;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}

			if((x + i) < colSize){
				int cell = lowy * colSize + x + i;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
		}
		if(upy < rowSize)
		{
			if((x - i) >= 0){
				int cell = upy * colSize + x - i;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}

			if((x + i) < colSize){
				int cell = upy * colSize + x + i;
				if(cellContain[cell] > 0)
					result->push_back(cell);
			}
		}
	}
}

void Grid::getObjectInCells(vector<int>* _cells, vector<int>* result)
{
	for(int i = 0; i < _cells->size(); i++)
	{
		if(nodeMap[_cells->at(i)].size() > 0){
			result->insert(result->end(),nodeMap[_cells->at(i)].begin(), nodeMap[_cells->at(i)].end());
		}
	}
}

// void Grid::getNodeInCells(vector<int>* _cells, vector<int>* result, int & objectNum)
// {
// 	for(int i = 0; i < _cells->size(); i++)
// 	{
// 		if(cellContain[_cells->at(i)] > 0){
// 			for(int j = 0; j < cellContain[_cells->at(i)]; j++){
// 				result->push_back(gridNodeMap[_cells->at(i)][j].nodeId);
// 				objectNum = objectNum + gridNodeMap[_cells->at(i)][j].objects.size();
// 			}
// 		}
// 	}
// }

inline int Grid::rightOf(int cid)
{
	if(cid  == -1)
		return -1;
	else
	{
		int result = cid + 1;
		return (result%colSize == 0) ? -1 : result;
	}
}

inline int Grid::leftOf(int cid)
{
	if(cid  == -1)
		return -1;
	else
	{
		int result = cid - 1;
		return (cid%colSize == 0) ? -1 : result;
	}
}

inline int Grid::upperOf(int cid)
{
	if(cid  == -1)
		return -1;
	else
	{
		int result = cid + colSize;
		return (result >= rowSize*colSize) ? -1 : result;
	}
}

inline int Grid::belowOf(int cid)
{
	if(cid  == -1)
		return -1;
	else
	{
		int result = cid - colSize;
		return (result<0) ? -1 : result;
	}
}

Grid::~Grid()
{
}


