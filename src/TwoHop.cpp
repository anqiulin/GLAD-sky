#include<cstdio>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include <sys/time.h>
#include<vector>
#include <xmmintrin.h>
#include<cmath>
#include"TwoHop.h"

TwoHop::TwoHop()
{
}

TwoHop::~TwoHop()
{
}

double TwoHop::GetTime()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

int TwoHop::LCAQuery(int _p, int _q){
	int p = toRMQ[_p], q = toRMQ[_q];

	if (p > q){
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;

	int i = LOGD[len], k = LOG2[len];

	q = q - i + 1;
	if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
		return RMQIndex[k][p];
	else return RMQIndex[k][q]; 
}

int TwoHop::distanceQuery(int p, int q){
	// added by Doris since the index start from 1 in hop index but start from 0 in graph index
	p = p + 1;
	q = q + 1;
	// end added
	if (p == q) return 0;
	int x = belong[p], y = belong[q];
	int result = 1;
	int lca = LCAQuery(x, y);
	if (lca == x || lca == y){
		if (lca == y){
			int v = y;
			y = x;
			x = v;
			v = p;
			p = q;
			q = v;
		}
		result =  dis[y][pos[x][posSize[x] - 1]];
	}
	else {

	
		int res = infinity;
		int *dx = dis[x], *dy = dis[y],*p2 = pos2[lca];
	//	_mm_prefetch(dx, _MM_HINT_T0);
	//	_mm_prefetch(dy, _MM_HINT_T0);
	//	_mm_prefetch(p2, _MM_HINT_T0);
		int ps = pos2Size[lca];
		for (int i = 0; i < ps; i++){
			int tmp = dx[p2[i]] + dy[p2[i]];
			if (res > tmp)
				res = tmp;
		}
		result = res;
	
	}
	
	return result;
	
	//return 0;
}

void TwoHop::scanIntArray(int *a, int n){
	fread(a, SIZEOFINT, n, fin);
}

int* TwoHop::scanIntVector(int *a){
	int _n;
	fread(&_n, SIZEOFINT, 1, fin);
	a = (int*)malloc(sizeof(int) * _n);
	scanIntArray(a, _n);
	return a;
}

void TwoHop::readIndex(char *file){

//	printf("*******************");

	double _time = GetTime();
	int tree_height = 0, tree_width = 0, most_sp = 0;
	fin = fopen(file, "rb");
	fread(&n, SIZEOFINT, 1, fin);
	int ts;
	fread(&ts, SIZEOFINT, 1, fin);
	TreeSize = ts;
	height = (int*)malloc(sizeof(int) * (ts + 1));
	for (int i = 0; i < ts; i++){
		fread(&height[i], SIZEOFINT, 1, fin);
	}
	belong = (int*)malloc(sizeof(int) * (n + 1));
	fread(belong, SIZEOFINT, n + 1, fin);
	toRMQ = (int*)malloc(sizeof(int) * (n + 1));
	fread(toRMQ, SIZEOFINT, n + 1, fin);
	int ris;
	fread(&ris, SIZEOFINT, 1, fin);
	fread(&ts, SIZEOFINT, 1, fin);
	EulerSeq = (int*)malloc(sizeof(int) * (ts + 1));
	RMQIndex = (int**)malloc(sizeof(int*) * (ris + 1));
	for (int i = 0; i < ris; i++){
		RMQIndex[i] = scanIntVector(RMQIndex[i]);
	}
	fread(&root, SIZEOFINT, 1, fin);
	cout << "root: " << root << endl;

	posSize = (int*)malloc(sizeof(int) * (n + 1));
	pos2Size = (int*)malloc(sizeof(int) * (n + 1));
	pos = (int**)malloc(sizeof(int*) * (TreeSize));
	pos2 = (int**)malloc(sizeof(int*) * (TreeSize));
	dis = (int**)malloc(sizeof(int*) * (TreeSize));
	chSize = (int*)malloc(sizeof(int) * (TreeSize));
	ch = (int**)malloc(sizeof(int*) * (TreeSize));

	for (int i = 0; i < TreeSize; i++){
		fread(&chSize[i], SIZEOFINT, 1, fin);
		ch[i] = (int*)malloc(sizeof(int) * chSize[i]);
		for (int j = 0; j < chSize[i]; j++){
			int x;
			fread(&x, SIZEOFINT, 1, fin);
			ch[i][j] = x;
		}
	}
	for (int i = 0; i < TreeSize; i++){
		int x;
		fread(&x, SIZEOFINT, 1, fin);
		fread(&posSize[x], SIZEOFINT, 1, fin);
		pos[x] = (int*)malloc(sizeof(int) * (posSize[x] + 1));
		fread(pos[x], SIZEOFINT, posSize[x], fin);
		if (posSize[x] > tree_width)
			tree_width = posSize[x];
		int _n;
		fread(&_n, SIZEOFINT, 1, fin);
		dis[x] = (int*)malloc(sizeof(int) * _n);
		fread(dis[x], SIZEOFINT, _n, fin);
		if (_n > tree_height)
			tree_height = _n;
	}
	printf("dis read finished!\n");
	for (int i = 0; i < TreeSize; i++){
		int x;
		fread(&x, SIZEOFINT, 1, fin);
		fread(&pos2Size[x], SIZEOFINT, 1, fin);
		pos2[x] = (int*)malloc(sizeof(int) * (pos2Size[x] + 1));
		fread(pos2[x], SIZEOFINT, pos2Size[x], fin);
		if (pos2Size[x] > most_sp)
			most_sp = pos2Size[x];
	}

	fclose(fin);
	printf("Load Index Time : %lf sec\n", (GetTime() - _time));
	printf("tree height: %d\n", tree_height);
	printf("tree width: %d\n", tree_width);
	printf("most search space: %d\n", most_sp);
}

void TwoHop::initLabel()
{
	LOG2 = (int*)malloc(sizeof(int) * (n * 2 + 10));
	LOGD = (int*)malloc(sizeof(int) * (n * 2 + 10));
	int k = 0, j = 1;
	for (int i = 0; i < n * 2 + 10; i++){
		if (i > j * 2){
			j *= 2;
			k++;
		}
		LOG2[i] = k;
		LOGD[i] = j;
	}
}

int TwoHop::H2HPath(int ID1, int ID2, vector<int>& vPath)//, vector<int>& vPathEdge)
{

	int p = ID1+1, q = ID2+1;
	int res = distanceQuery(p-1, q-1);
	vPath.push_back(p-1);
	int tn = p;
        while (p != q){ 
        	int pq = distanceQuery(p-1, q-1); 
                for (int i = 0; i < Degree[p]; i++){
  			int x = Neighbor[p][i];
                        int xq = distanceQuery(x-1, q-1); 
                        if (xq + Weight[p][i] == pq){
 	                       p = x;
                               vPath.push_back(p-1);
		       	       //vPathEdge.push_back(adjListEdge[tn-1][i].second);
			    // cout<<tn-1<<": "<<p-1<<" = "<<adjListEdge[tn-1][i].second<<endl;
			       tn = x; 
                               break;
                         }
                 }
        }
        return res;
}

void TwoHop::readgraph(char *filename){
		FILE * file = fopen(filename, "r");
		int n, m;
		//cout << filename <<endl;
		//fstream c("testData/graph.txt",ios::out);
		//c <<"111";
		//cout << file <<endl;
		fscanf(file, "%d %d", &n, &m);
		Degree = (int*)malloc(sizeof(int) * (n + 1));
		vector< vector<pair<int, int> > > nb;
		vector<pair<int, int> > v;
		v.clear();
		for (int i = 0; i <= n; i++){
		//	Degree[i] = 0;
			nb.push_back(v);	
		}
	//	cout << n << " " << m << endl;
		for (int i = 0; i < m; i++){
			int x, y, z;
			fscanf(file, "%d %d %d", &x, &y, &z);
	//		Degree[x]++;
	//		cout << x << " " << y << " " << z << endl;
			nb[x].push_back(make_pair(y, z));
		}
		Neighbor = (int**)malloc(sizeof(int*) * (n + 1));
		Weight = (int**)malloc(sizeof(int*) * (n + 1));
		for (int i = 1; i <= n; i++){
			Degree[i] = nb[i].size();
			Neighbor[i] = (int*)malloc(sizeof(int) * nb[i].size());
			Weight[i] = (int*)malloc(sizeof(int) * nb[i].size());
			for (int j = 0; j < nb[i].size(); j++){
				Neighbor[i][j] = nb[i][j].first;
				Weight[i][j] = nb[i][j].second;
			}
		}
	}