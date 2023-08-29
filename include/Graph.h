/*
	Author: Qiulin An
	Date: 6/5/2023
*/
#ifndef GRAPH_H
#define GRAPH_H
#include<vector>
#include<map>
#include <mutex> 
#include <atomic>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Point.h"
#include<bitset>
using namespace std;

 //NY-Density: 2%-158 4%-317 8%-634	16%-1268 32%-2537
// NY-KIND: 3%-634 6%-1268 9%-1903 12%-2537 15%-3172
static const int KIND = 158;
/*
	This graph class suports incrementally build and dynamic add and delete objects on the nodes.
	However the dynamic construction of the graph is not available.

	Modified the graph such that the id of a node corresponds to the index of the node in the vector
	Before add nodes and edges of the graph, we need to sort the nodes and edges in advance.
	After sorting, nodes are added one by one, and all the outgoing edges of a node are added right after
	the node is added.
*/

class Graph
{
public:
	bool isDirected;
	int numNode;
	int numEdge;
	int numObject;
	struct Vertex {
		int objectCount = 0;
		vector<int> objectList;
	};
	vector<Vertex> nodes;
	/*
	//set KEYS !!!!!!!
	bitset<KIND>Query;
	vector<vector<int>> RSP;// each keyword contain pois (0~299)
	vector<map<int, int>>KEYS;//store each poi's keywords using map, and node start from zero
	void set_nodeKEYS_and_Query(string filename, vector<int>& querykey);//init Query and each node keywors using map
	*/
	vector<bitset<KIND> > NodesBit;
	vector<int> POIObjects;
	bitset<KIND> Query;

	int pointer;
	vector<int> nodeIndex; // which is a numNode * 2 array
	vector<int> neighbor; // the size equals to the numEdge * 2
	vector<int> movingObjects;
	Graph();
	~Graph();
	Graph(Graph* graph);
	void clear();
	void init(int _numNode, int _numEdge, int _pointer, int _numObeject);
	void addEdge(int _node1, int _node2, int _weight);
	void addObject(int _node, int _object);
	void deleteObject(int _object);
	void updateObject(int _nodeNew, int _object);
	void buildGraph(vector<int*> _edges, int _numberNode, vector<int> _movingObejcts);
	vector<int*> readGraph(string _filePath);
	vector<Point> readNodeCoordinate(string _filePath, int _numPoint);
	vector<int*> readGraphBJ(string _filePath);
	vector<Point> readNodeCoordinateBJ(string _filePath, int _numPoint);
	vector<int> readMovingObjects(string _filePath);
	vector<int> readQuery(string _filePath);
};
#endif
