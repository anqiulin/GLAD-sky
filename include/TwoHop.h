#ifndef TWOHOP_H
#define TWOHOP_H
#include <algorithm>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include<sys/time.h>
#include<vector>
#include<xmmintrin.h>
#include<cmath>

using namespace std;
class TwoHop{

public:

	const int infinity = 999999999;
	const int SIZEOFINT = 4;

	int *toRMQ, *height, **RMQIndex;
	int *belong;
	int root, TreeSize;
	int **rootToRoot, *rootSite;
	int **dis, **pos, **pos2;
	int *posSize, *pos2Size;
	int *chSize;
	int ** ch;
	int *LOG2, *LOGD; 
	int rootSize;
	int *Degree;
	int **Neighbor, **Weight;

	int counter = 0;
	double RunTime = 0;

	FILE *fin;
	int n;
	int *EulerSeq;

	TwoHop();
	~TwoHop();

	double GetTime();
	void initLabel();
	void readIndex(char *file);
	int* scanIntVector(int* a);
	void scanIntArray(int* a, int n);
	int distanceQuery(int p, int q);
	int LCAQuery(int p, int q);
	int H2HPath(int ID1, int ID2, vector<int>& vPath);
	void readgraph(char *filename);
};
#endif
