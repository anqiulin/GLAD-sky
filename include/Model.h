#include<vector>
#include<mutex>
#include<atomic>
#include<algorithm>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<pthread.h>
#include"Graph.h"
#include"Grid.h"
#include"MyHeap.h"
#include"TwoHop.h"
#include"GlobalVar.h"
#include<time.h>
#include<unordered_map>
#include<math.h>
#include<cmath>
 #include <limits.h>
using namespace std;



vector<int> GreedyAglorByDis(TwoHop* _twoHop, Graph* _graph, vector<int>node, int s, bitset<KIND> qu);
vector<vector<int>> getcombine(int n, int m);
void fun(int level);
void func(int n, int m);
int GreedyAglorByVS(Graph* _graph, vector<int>node, bitset<KIND> query);

double getTime(){
	double time;
	struct timeval tp;
	gettimeofday(&tp, NULL);
	int tmp =(long long)tp.tv_sec % 86400;
	time = tmp * 1000000 + (double)tp.tv_usec;
	return time;
}



bool thisCompare(INT2 a, INT2 b){return (a.elem[0] < b.elem[0]);}




void readPOI(string _filePath, vector<int>* _poi) {
	ifstream poiFile(_filePath);
	if (!poiFile) {
		printf("readPOI::Warning: cannot open the file! \n");
		return;
	}
	string line;
	while (poiFile && getline(poiFile, line)) {
		_poi->push_back(stoi(line, nullptr, 0) - 1);
	}

}



void initGLADkNN(Graph* graph, vector<Point>* pointList, Grid* grid, TwoHop* twoHop)
{
	printf("Start the Program! \n");
	vector<int*> edgesList;
	if(DATA == 0)
		edgesList = graph->readGraph("./map/data_NY/tol/NY.edge");
	else if(DATA == 1)
		edgesList = graph->readGraph("./map/data_cal/tol/cal.edge");
	else if(DATA == 2)
		edgesList = graph->readGraph("./map/data_COL/tol/COL.edge");

	printf("Finish the loading of graph! \n");
	int edgeSize = edgesList.size();
	printf("edge.size = %d \n", edgeSize);
	
	// vector<int> movingObjects(NUMOBJECT);
	// generateQueryListUniform(NUMOBJECT, &movingObjects);
	// printf("movingObject.size = %d \n", (int)movingObjects.size());
	// ofstream out;
	// string outfile = "./movingObject-100000";
	// out.open(outfile, ios_base::out);
	// for(int i = 0; i < NUMOBJECT; i++)
	// 	out << i << "," <<  movingObjects[i] << "\n";

	MBR* extend = new MBR();
	for(int i = 0; i < NUMPOINT; i++)
		extend->unionWith(&pointList->at(i));
	Point** objectPointList = new Point*[graph->POIObjects.size()];
	for(int i = 0; i < graph->POIObjects.size(); i++){
		objectPointList[i] = &pointList->at(graph->POIObjects[i]);
		objectPointList[i]->id = i;
		objectPointList[i]->nodeId = graph->POIObjects[i];
	}

	grid->extend = extend;
	grid->build(objectPointList, graph, CELLSIZE);
	cout<<"cellsize = "<<NUMOBJECT<<" "<<CELLSIZE<<endl;
	printf("Finish the building of grid! \n");
	graph->buildGraph(edgesList, NUMPOINT, graph->POIObjects);
	edgesList.clear();
	for(int i = 0; i < graph->POIObjects.size(); i++) {
		ObjectInfo* tmp = new ObjectInfo;
		tmp->availability = 1;
		OBJECTLIST.push_back(tmp);
	}
	graph->POIObjects.clear();
	

	cout<<"ready index"<<endl;
	string filepath_idex, filepath_edge;
	if(DATA == 0)
	{
		filepath_idex = "./map/data_NY/h2h/NY.index";
		filepath_edge = "./map/data_NY/h2h/NY.edge";
	}
	else if(DATA == 1)
	{
		filepath_idex = "./map/data_cal/h2h/cal.index";
		filepath_edge = "./map/data_cal/h2h/cal.edge";
	}
	else if(DATA == 2)
	{
		filepath_idex = "./map/data_COL/h2h/COL.index";
		filepath_edge = "./map/data_COL/h2h/COL.edge";
	}
	char *file_idex = new char[filepath_idex.length() + 1];
	char *file_edge = new char[filepath_edge.length() + 1];
	strcpy(file_idex, filepath_idex.c_str());
	strcpy(file_edge, filepath_edge.c_str());
	twoHop->readIndex(file_idex);
	twoHop->readgraph(file_edge);
	twoHop->initLabel();
	printf("Finish the construction of graph! \n");
}

void getPOIcandidate(Graph* _graph, Grid* _grid, vector<Point>* _pointList, int s, int t, bitset <KIND> query, vector<int>& PoiCand)
{
	vector<int> cell;
	// bitset<KIND> cellBit;
	cell.clear();
	
	// _grid->getCell_MBR(&_pointList->at(s), &_pointList->at(t), query, &cell);
	
	// for(auto c : cell)
	// {
	// 	cout << c << " ";
	// }
	// cout << endl;

	// cell.clear();
	// cout << "------------";
	_grid->getcell_Straight(&_pointList->at(s), &_pointList->at(t), query, &cell);
	
	// for(auto c : cell)
	// {
	// 	cout << c << " ";
	// }
	// cout << endl;
	vector<int> tmp;
	tmp.clear();

	_grid->getPoiInCells(&cell, &tmp);

	for(int i = 0; i < tmp.size(); i++) //delete nodes without querykey
	{
		if((_graph->NodesBit[tmp[i]] & query).count() != 0)
		{
			PoiCand.push_back(tmp[i]);
		}
	}
}

vector<vector<int>>a;
vector<int>len;
vector<int>p;
int count_combine=0;

vector<vector<int>> res_combine;
vector<int> tmp_combine;

vector<vector<int>> expand_enum(TwoHop* _twohop, Graph* _graph, vector<int> node, vector <int> qu, bitset <KIND> query, int s, int t)
{
	vector<vector<int>> KeyNum; //Keywords number index
	vector<vector<int>> RSPOI; //Keywords inverted index
	RSPOI.clear();
	KeyNum.clear();

	

	for(int i = 0; i <= qu.size(); i++)
	{
		vector<int> tmp;
		// tmp.clear();
		KeyNum.push_back(tmp);
	}
	for(int i = 0; i < qu.size(); i++)
	{
		vector<int> tmp;
		RSPOI.push_back(tmp);
	}
	for(int i = 0; i < node.size(); i++)
	{
		KeyNum[(_graph->NodesBit[node[i]] & query).count() ].push_back(node[i]);
		for(int j = 0; j < qu.size(); j++)
		{
			if((_graph->NodesBit[node[i]] & query).count() == 1) 
			{
				if(_graph->NodesBit[node[i]].test(qu[j]))
				{
					RSPOI[j].push_back(node[i]);
				}
			}
			
		}
	}
	// cout << "RSPOI[0][0]" << endl;
	// for(auto c : RSPOI)
	// {
	// 	for(auto k : c)
	// 	{
	// 		cout << k << " " ;
	// 	}
	// 	cout << endl;
	// }
	// cout << "key" << endl;
	// int pp = 0;
	// for(auto c : KeyNum)
	// {
	// 	cout << pp++ << endl;
	// 	for(auto k : c)
	// 	{
	// 		cout << k << " " ;
	// 	}
	// 	cout << endl;
	// }
	

	vector<vector<int> > ret;

	vector<int>temp;
	temp = GreedyAglorByDis(_twohop,_graph,node,s,query);

	int vs_lower = GreedyAglorByVS(_graph,node,query);

	ret.push_back(temp);
	int vs_upper = temp.size();
	if(vs_lower == vs_upper)
		return ret;
	for(int i = vs_upper-1; i >= vs_lower; i--)
	{
		vector<vector<int>> combine = getcombine(i,qu.size());

		for(auto tmp : combine)
		{
			int flag = 0;
			for(auto c : tmp)
			{
				if(KeyNum[c].size() == 0)
				flag = 1;
				// cout << c << " ";
			}
			// cout << endl;
			if(flag)
				continue;
			a.clear();
			len.clear();
			p.clear();
			res_combine.clear();

			int Cnum;
			for(Cnum = 0; Cnum < tmp.size(); Cnum++)
			{
				if(tmp[Cnum]!=1)
				{
					a.push_back(KeyNum[tmp[Cnum]]);
					len.push_back(KeyNum[tmp[Cnum]].size());
					p.push_back(0);
				}
				else
					break;
			}

			count_combine = Cnum;

			vector<vector<int>> tmpres;

			if(count_combine == 1)
			{
				for(int j = 0; j < a[0].size(); j++)
				{
					tmpres.push_back({a[0][j]});
				}
			}
			else{
				fun(0);
			
				for(auto c : res_combine)
				{
					temp.clear();
					temp.insert(temp.begin(),c.begin(),c.end());
					tmpres.push_back(temp);
				}
			}
			
			for(auto c : tmpres)
			{
				a.clear();
				len.clear();
				p.clear();
				res_combine.clear();

				bitset<KIND>t;
				t.reset();
				for(auto k : c)
				{
					a.push_back({k});
					len.push_back(1);
					p.push_back(0);
					t |= _graph->NodesBit[k] & query;
				}
				if((t.count() + tmp.size() - count_combine) != qu.size())
					continue;

				for(int j = count_combine; j < tmp.size(); j++)
				{
					for(int k = 0; k < RSPOI.size(); k++)
					{
						// cout << RSPOI[k][0] << endl;
						if(RSPOI[k].size() == 0)
							continue;
						bitset<KIND> cnt(_graph->NodesBit[RSPOI[k][0]] & query);

						if((cnt & (~t)).count() > 0)
						{
							t |= cnt;
							a.push_back(RSPOI[k]);
							len.push_back(RSPOI[k].size());
							p.push_back(0);
							break;
							
						}
					}
				}

				if(t.count() != query.count())
					continue;

				count_combine = i;
				
				fun(0);
				
				for(auto k : res_combine)
				{
					bitset<KIND>cnt;
					cnt.reset();
					temp.clear();
					for(auto m : k)
					{
						cnt |= _graph->NodesBit[m] & query;
						temp.push_back(m);
					}
					if(cnt.count() == query.count())
						ret.push_back(temp);
				}
				                                                                                  
			}
		}
		
	}
	cout << "ret size" << ret.size() << endl;
	// for(auto c : ret)
	// {
	// 	for(auto k : c)
	// 	{
	// 		cout << k << " " ;
	// 	}
	// 	cout << endl;
	// }
	return ret;
}
vector<vector<int>> expand_wpoi(TwoHop* _twohop, Graph* _graph, vector<int> node, vector <int> qu, bitset <KIND> query, int s, int t)
{

	vector<vector<int>> RSPOI;
	RSPOI.clear();

	for(int i = 0; i < qu.size(); i++)
	{
		vector<int> tmp;
		// tmp.clear();
		RSPOI.push_back(tmp);
	}

	for(int i = 0; i < node.size(); i++)
	{
		for(int j = 0; j < qu.size(); j++)
		{ 
			if(_graph->NodesBit[node[i]].test(qu[j]))
			{
				RSPOI[j].push_back(node[i]);
			}
		}
	}

	unordered_map<int,pair<double,double>> p;
	p.clear();
	for(auto &k : node)
	{
		double pdist =  (double)(((_twohop->distanceQuery(s, k) + _twohop->distanceQuery(k, t)) - _twohop->distanceQuery(s, t)) *1.0 / _twohop->distanceQuery(s, t) *1.0);
		bitset <KIND> tmp(query & _graph->NodesBit[k]);

		double pstop = (double)(((query.count() - tmp.count()) *1.0) / (query.count() *1.0));
		p[k] = make_pair(pdist, pstop);
	}
	
	// fstream out("plist.txt",ios::app);
	// out << "s:" << s << " t:" << t << endl;
	// for(auto c : p)
	// {
	// 	out << "id:" << c.first <<  " dis:" << c.second.first << " pstop:" << c.second.second << endl;  
	// }

	vector<double> w;
	w.clear();

	for(int i = 0; i < 1000; i++)
	{
		double t = double (i*1.0 /1000.0);
		w.push_back(t);
	}

	vector<vector<int>> ret;
	ret.clear();

	for(int wp = 0;  wp < w.size(); wp++)
	{ 
		vector <int> tmp;
		tmp.clear();
		bitset <KIND> t;
		t.reset();
		for(int i = 0; i < RSPOI.size(); i++)
		{
			if(t.test(qu[i]))
				continue;
			double min_p= 100.0;
			int id = 0;
			for(int j = 0; j < RSPOI[i].size(); j++)
			{
				double pp = p[RSPOI[i][j]].first * w[wp] + p[RSPOI[i][j]].second * (1.0-w[wp]);
				if(pp < min_p)
				{
					min_p = pp;
					id = RSPOI[i][j];
				}
			}
			tmp.push_back(id);
			t |= _graph->NodesBit[id] & query;
			if(t.count() == query.count())
				break;
		}
		ret.push_back(tmp);
	}
	// out << "------------" << endl;
	// for(auto c : ret)
	// {
	// 	for(auto k : c)
	// 	{
	// 		out << k << " ";
	// 	}
	// 	out << endl;
	// }

	return ret;

}

int getTSP(TwoHop* _twohop, int s,int t,vector<int> node, vector<int> &order){
	int dist=0;
	int index = 0;
	order.clear();
	order.push_back(s);
	while(node.size()> 0){
	int min=INT_MAX;
		for(int i=0;i<node.size();i++){
			if(min > _twohop -> distanceQuery(s,node[i])){
				min = _twohop -> distanceQuery(s,node[i]);		
				index = node[i];
			}
		}
		dist+=min;
		order.push_back(index);
		for(vector<int>::iterator iter=node.begin();iter!=node.end();iter++){        //从vector中删除指定的某一个元素 
    			if(*iter==index){
        			node.erase(iter);
        			break;
    			}
		}
		s = index;	
	}
	dist+= _twohop -> distanceQuery(index,t);
	order.push_back(t);

	return dist;
}

void getFullPath(TwoHop* _twohop, vector<int> order, vector<int> &path)
{
	path.push_back(order[0]);
	for(int i = 0; i < order.size()-1; i++)
	{
		vector<int> tmp;
		tmp.clear();
		int s = _twohop->H2HPath(order[i],order[i+1],tmp);
		path.pop_back();
		path.insert(path.end(),tmp.begin(),tmp.end());
	}
}

bitset<KIND> set_str_to_binary(string ss){
	bitset<KIND>binary;
		char str[200]; //according to the keywordstr.length
			strcpy(str,ss.c_str());
			const char * split = ",";
    		char * p;
    		p = strtok (str,split);
    		while(p!=NULL) {
			        if(strcmp(p,"-1") ==0)
						break;
					binary.set(atoi(p));   
        		p = strtok(NULL,split);
      		}
        return binary;
}

vector<int> set_str_to_vec(string ss){

	vector<int> ret;
	ret.clear();
		char str[100]; //according to the keywordstr.length
			strcpy(str,ss.c_str());

			const char * split = ",";
    		char * p;
    		p = strtok (str,split);

    		while(p!=NULL) {
			        if(strcmp(p,"-1") ==0)
						break;
					ret.push_back(atoi(p));   
        		p = strtok(NULL,split);
        }
        return ret;
}

void GLADsearch(Graph* _graph, Grid* _grid, TwoHop* _twohop, vector<Point>* _pointList, int s, int t, string querystr, vector<PA>& Res, vector<vector<int>>& path)
{
	// ofstream out("POICandiate1.txt",ios::app);
	Res.clear();

	bitset<KIND> query = set_str_to_binary(querystr);
	vector<int> qu = set_str_to_vec(querystr);


	vector<int>Poicand;
	Poicand.clear();
	getPOIcandidate(_graph, _grid, _pointList, s, t, query, Poicand);

	// out << "s: " << s << " t:" << t << endl;
	// for(auto c : Poicand)
	// {
	// 	out << c << " ";
	// }
	// out << endl;

	vector<vector<int>> candidate;
	candidate.clear();

	// candidate = expand_enum(_twohop, _graph, Poicand, qu, query, s, t);
// 
	candidate = expand_wpoi(_twohop, _graph, Poicand, qu, query, s, t);

	
	for(int i = 0; i < candidate.size(); i++)
	{
		vector<int> route;
		route.clear();

		// std::chrono::high_resolution_clock::time_point t3;
		// std::chrono::high_resolution_clock::time_point t4;
		// std::chrono::duration<double> time_span;
		// t3 = std::chrono::high_resolution_clock::now();	
		int dis = getTSP(_twohop, s, t, candidate[i], route);
		// t4 = std::chrono::high_resolution_clock::now();	
		// time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t4 - t3);
		// cout << "gettsp_time: " << time_span.count() << endl;

		if(Res.empty())
		{
			Res.push_back(PA(route, dis, route.size()-2));
			continue;
		}
		
		while ((Res.back().dis > dis && Res.back().park >= (route.size()-2)) || (Res.back().dis >= dis && Res.back().park > (route.size()-2)))
		{
			Res.pop_back();
			if(Res.empty())
				break;
		}
		
		if(Res.empty() || (Res.back().dis < dis && Res.back().park > (route.size()-2)) || (Res.back().dis > dis && Res.back().park < (route.size()-2)))
			Res.push_back(PA(route, dis, route.size()-2));
	}
	for(int  i = 0; i < Res.size(); i++)
	{
		vector<int> tmp;
		tmp.clear();
		getFullPath(_twohop, Res[i].pathnode, tmp);
		path.push_back(tmp);
	}
}

vector<int> GreedyAglorByDis(TwoHop* _twoHop, Graph* _graph, vector<int>node, int s, bitset<KIND> qu){

	int maxKey = 0;
	bitset<KIND>q;

	bitset<KIND>temp;
	vector<int>nodeSet;
	q|=qu;

	int len = INT_MAX;
	int minlen = 0;
	while(q.count()>0){
		maxKey = -1;
		for(int i=0;i<node.size();i++){
			temp.reset();
			temp|=(q&_graph->NodesBit[node[i]]);//取共同 
			if(temp.count()>0 || maxKey==-1){
				if(maxKey==-1){
					maxKey = node[i];
					minlen = _twoHop->distanceQuery(s,node[i]);
					// outfile_test << "2" << endl;
					continue;
				}
				len = _twoHop->distanceQuery(s,node[i]);
					//cout<<node[i]<<" no "<<len<<endl;
					if(len<minlen){
						maxKey = node[i];
						minlen = len;
						// outfile_test << "3" << endl;
					//cout<<node[i]<<" ye "<<minlen<<endl;
					}
				
			}
		}
		if(maxKey!=-1){
			temp.reset();
			temp|=(q&_graph->NodesBit[maxKey]);//取共同
			if(temp.count() > 0)
			{
				nodeSet.push_back(maxKey);
				// outfile_test << "4" << " " << maxKey << endl;
				s = maxKey;
				//cout<<"////"<<maxKey<<endl;
				bitset<KIND>t;
				//bitset<20>t;
				t|=_graph->NodesBit[maxKey]; 
				q&=(~t);
				for(vector<int>::iterator iter=node.begin();iter!=node.end();iter++){  //从vector中删除指定的某一个元素 
    				if(*iter==maxKey){
        				node.erase(iter);
        				break;
    				}
				}
			}
			else{  //aql modifiy
				for(vector<int>::iterator iter=node.begin();iter!=node.end();iter++){  //从vector中删除指定的某一个元素 
    				if(*iter==maxKey){
        				node.erase(iter);
        				break;
    				}
				}
			}	
		}
		// outfile_test << "q:" << q.to_string() << endl;

	}
	return nodeSet;
}

int GreedyAglorByVS(Graph* _graph, vector<int>node, bitset<KIND> query)
{
	bitset<KIND> cnt(query);
	vector<int> Minpath;
	bitset<KIND> temp;
	Minpath.clear();
	int maxKey;
	int maxnum;
	int Num = 0;
	while(cnt.count() > 0)
	{
		maxKey = -1;
		maxnum = -1;
		for(int i = 0; i < node.size(); i++)
		{
			temp.reset();
			temp |= (cnt & _graph->NodesBit[node[i]]);
			if(temp.count() > 0)
			{
				Num = temp.count();
				if(Num > maxnum){
					maxnum = temp.count();
					maxKey = i;
				}
			}
		}
		if(maxKey!=-1){
			temp.reset();
			temp |= (cnt & _graph->NodesBit[node[maxKey]]);//取共同
			if(temp.count() > 0)
			{
				Minpath.push_back(node[maxKey]);
				bitset<KIND>t;
				t |= _graph->NodesBit[node[maxKey]]; 
				cnt &= (~t);
				node.erase(node.begin()+maxKey);
			}
		}
	}
	return Minpath.size();
}

vector<int>record;
vector<vector<int>> result;

void func(int n, int m)
{
    if (n == 2)
    {
        for (int i = 1; i < m; ++i)
        {
            record.push_back(i);
            record.push_back(m - i);
            result.push_back(record);
            record.pop_back();
            record.pop_back();
        }

        return;
    }
    for (int beg = 1; beg <= m - 2; ++beg)
    {
        record.push_back(beg);
        func(n - 1, m - beg);
        record.pop_back();
    }

}

vector<vector<int>> getcombine(int n, int m)
{
	record.clear();
	result.clear();
	func(n, m);
    vector <vector<int>>res;
    for (auto tmp : result)
    {
        sort(tmp.begin(), tmp.end(), greater<int>());
        if (res.empty())
        {
            res.push_back(tmp);
            continue;
        }
        int flag = 1;
        for (int j = 0; j < res.size(); j++)
        {
            if (tmp == res[j])
            {
                flag = 0;
                break;
            }
        }
        if (flag)
            res.push_back(tmp);
    }
    return res;
}

void fun(int level) {
    int i;
    if (level < count_combine) {
        for (i = 0; i < len[level]; i++) {
            p[level] = i;
            fun(level + 1);
        }
    }
    else {
        tmp_combine.clear();
        for (i = 0; i < count_combine; i++) {
            tmp_combine.push_back(a[i][p[i]]);
        }
        res_combine.push_back(tmp_combine);
    }
}