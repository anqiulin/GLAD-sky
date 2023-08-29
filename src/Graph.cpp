/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/

#include "Graph.h"
vector<string> split(string s, char delim);

Graph::Graph()
{
}

Graph::~Graph()
{
}

Graph::Graph(Graph* graph)
{
	init(graph->numNode, graph->numEdge, graph->pointer, graph->numObject);
	nodes = graph->nodes;
	nodeIndex = graph->nodeIndex;
	neighbor = graph->neighbor;
	movingObjects = graph->movingObjects;
}

void Graph::clear() {
	numNode = numEdge = pointer = 0;
	neighbor.clear();
	nodeIndex.clear();
	nodes.clear();
	movingObjects.clear();
}

void Graph::init(int _numNode, int _numEdge, int _pointer, int _numObject)
{
	clear();
	isDirected = true;
	numNode = _numNode;
	numEdge = _numEdge;
	pointer = _pointer;
	numObject = _numObject;
	nodeIndex.assign(_numNode * 2, -1);
	neighbor.assign(_numEdge * 2, -1);
	movingObjects.assign(_numObject, -1);
	Vertex elem;
	nodes.assign(_numNode, elem);
}

void Graph::addEdge(int _node1, int _node2, int _weight)
{
	if (nodeIndex[_node1 * 2] == -1) {
		nodeIndex[_node1 * 2] = pointer;
		nodeIndex[_node1 * 2 + 1] = 1;
	}
	else {
		nodeIndex[_node1 * 2 + 1]++;
	}
	neighbor[pointer * 2] = _node2;
	neighbor[pointer * 2 + 1] = _weight;
	pointer++;
}


void Graph::buildGraph(vector<int*> _edges, int _numNode, vector<int> _movingObejcts) {
	init(_numNode, _edges.size(), 0, _movingObejcts.size());
	for (int i = 0; i < _edges.size(); i++) {
		addEdge(_edges[i][0], _edges[i][1], _edges[i][2]);
	}
}


vector<int*> Graph::readGraph(string _filePath) {
	vector<int*> _edgesList;
	ifstream graphFile(_filePath);
	if (!graphFile) {
		printf("readGraph::Warning: cannot open the file! \n");
		return _edgesList;
	}
	string line;
	while (graphFile && getline(graphFile, line)) {
		vector<string> rowdata = split(line, ' ');
		int* coordinates = new int[3];
		coordinates[0] = stoi(rowdata[0], nullptr, 0) - 1;
		coordinates[1] = stoi(rowdata[1], nullptr, 0) - 1;
		coordinates[2] = stoi(rowdata[2], nullptr, 0);
		_edgesList.push_back(coordinates);
		rowdata.clear();
	}
	return _edgesList;
}

vector<Point> Graph::readNodeCoordinate(string _filePath, int _numPoint) {
	for(int i=0;i<_numPoint+1;i++)
	{
		bitset<KIND>tmp;
		NodesBit.push_back(tmp);
		// NodesWord.push_back(-1);
	}
	POIObjects.clear();
	vector<Point> pointList;
	Point tmp;
	pointList.assign(_numPoint, tmp);
	ifstream nodeFile(_filePath);
	if (!nodeFile) {
		printf("readNodeCoordinate::Warning: cannot open the file! \n");
		return pointList;
	}
	string line;
	int pointId = 0;
	while (nodeFile && getline(nodeFile, line)) {
		vector<string> rowdata = split(line, ' ');
		pointList[pointId].id = pointId;
		pointList[pointId].coord[0] = stoi(rowdata[1], nullptr, 0);
		pointList[pointId].coord[1] = stoi(rowdata[2], nullptr, 0);
		vector<string> keys = split(rowdata[3], ',');
		if(stoi(keys[0]) == -1)
		{
			NodesBit[pointId].reset();
		}
		else
		{
			POIObjects.push_back(pointId);
			for(int i = 0; i < keys.size(); i++)
			{
				int key = stoi(keys[i]);
				if(key == -1)
					break;
				NodesBit[pointId].set(key);
			}
		}
		pointId++;
		keys.clear();
		rowdata.clear();
	}
	return pointList;
}


vector<string> split(string s, char delim) {
	stringstream ss(s);
	string item;
	vector<string> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}
/*
void Graph::set_nodeKEYS_and_Query(string filename, vector<int>& querykey){
	//set Query
	for(int k:querykey){
		Query.set(k);
	}
	cout<<"QueryNum = "<<Query.count()<<endl;
	//set RSP
	for(int i = 0;i < KIND; i++){
		vector< int >tmp;
		RSP.push_back(tmp);
	}

	//set nodeKEYS
	fstream fp(filename);
        int nid;
        string ss;
        long x,y; 
        while(fp >> nid >> x >> y >> ss){
        ///	cout<<nid<<" "<<ss<<endl;
                char str[500]; //according to the keywordstr.length (50*4)
                strcpy(str,ss.c_str());
                const char * split = ",";
                char * p;
                p = strtok (str,split);
                map<int,int>mp;
                mp.clear();
                while(p!=NULL) {
                                if(strcmp(p,"-1") ==0)
                                        break;
                                mp.insert(make_pair(atoi(p),1));
                                RSP[atoi(p)].push_back(nid-1);
                                p = strtok(NULL,split);
                }
                KEYS.push_back(mp);
        }
	cout<<"KEYS map and Query init complete!"<<endl;
        fp.close();

}
*/
