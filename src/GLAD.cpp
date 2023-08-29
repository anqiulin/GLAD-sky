/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cilk/cilk.h>
#include <sys/time.h>

#include"Graph.h"
#include"Grid.h"
#include"Point.h"
#include"MBR.h"
#include"Model.h"
#include"TwoHop.h"

using namespace std;


int main()
{
	Graph* graph = new Graph();
	Grid* grid = new Grid();
	TwoHop* twoHop = new TwoHop();
	vector<Point> pointList;

	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	
	std::chrono::duration<double> time_span;
	std::chrono::duration<double> GLAD_time_totol;

	if(DATA == 0)
	{
		NUMPOINT = 264346;
		pointList = graph->readNodeCoordinate("./map/data_NY/tol/Density2/NY_keys_15_3%.node", NUMPOINT);
	}
	else if(DATA == 1)
	{
		NUMPOINT = 21048;
		pointList = graph->readNodeCoordinate("./map/data_cal/tol/Density8/cal_keys_15_3%.node", NUMPOINT);
	}
	else if(DATA == 2){
		NUMPOINT = 435666;
		pointList = graph->readNodeCoordinate("./map/data_COL/tol/Density8/COL_keys_15_3%.node", NUMPOINT);
	}


	initGLADkNN(graph, &pointList, grid, twoHop);

	string querystr;

	// querystr= "28,33"; //cal
	// querystr= "10,16,23,32"; //4
	// querystr= "8,18,24,17,28,26"; //6
	// querystr= "7,10,11,24,26,31,38,45"; //8
	// querystr= "6,18,24,33,46,2,34,25,47,19"; //10

	// querystr= "156,343";//NY
	// querystr= "147,245,354,498";
	// querystr= "99,154,295,334,347,456";
	// querystr= "137,247,359,468,521,570,591,621";
	// querystr= "45,89,145,169,258,290,358,402,479,587";

	// querystr= "447,450"; //1000COL
	// querystr= "178,379,511,479"; //4
	// querystr= "278,384,708,930,528,875"; //6
	// querystr= "115,551,247,330,317,368,443,741"; //8
	// querystr= "978,349,407,146,421,690,643,726,829,626"; //10

	// querystr= "426,635,727,730,271,556"; //6%
	// querystr= "598,676,1047,1196,1410,1411"; //9%
	// querystr= "936,1265,1417,1589,1603,1769"; //12%
	// querystr= "453,809,1505,2057,2411,2585"; //15%

	querystr= "11,58,57,96,137,141"; //Density:2%
	// querystr= "11,95,127,216,225,272";
	// querystr= "287,369,443,576,825,975";
	// querystr= "362,428,946,1123,1411,1738";

	string Scout = "COL keywords:6%  querykey:6-{6,18,24,33,46,2,34,25,47,19}  Dis:0~20% Path generation: enum";

	fstream infile("testData/NY_tol/query/query_1000_0-20%.txt",ios::in);

	ofstream outfile("testData/NY_tol/analysis_res/Density2/3%_6_0-20%_GLAD_weight.txt",ios::out);
	ofstream outfile_res("testData/NY_tol/result/Density2/result.txt",ios::app);
	
	// fstream infile("testData/Cal_tol/query/query_1000_0-20%.txt",ios::in);

	// ofstream outfile("testData/Cal_tol/analysis_res/Density8/3%_6_0-20%_GLAD_weight.txt",ios::out);
	// ofstream outfile_res("testData/Cal_tol/result/Density8/result.txt",ios::app);

	// fstream infile("testData/COL_tol/query/query_1000_0-20%.txt",ios::in);

	// ofstream outfile("testData/COL_tol/analysis_res/Density8/3%_10_0-20%_GLAD_enum1.txt",ios::out);
	// ofstream outfile_res("testData/COL_tol/result/Density8/result.txt",ios::app);


	int s, t, len;
	int num = 0;
	srand (time(NULL));

	while(infile >> s >> t >> len)
	{
		// if(num < 73)
		// 	{
		// 		num++;
		// 		continue;
		// 	}		
		if(num == 1000)
			break;
		num++;
		cout << "num:" << num << " s:" << s << " t:" << t << endl; 
		vector<PA> Res;
		Res.clear();
		// int s=154759, t=211728;
		vector<vector<int>> path;
		path.clear();

		t1 = std::chrono::high_resolution_clock::now();
		GLADsearch(graph, grid, twoHop, &pointList, s, t, querystr, Res, path);
		t2 = std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
		GLAD_time_totol += time_span;
		cout << "base times: " << time_span.count() << endl; 
		cout << "totol times: " <<GLAD_time_totol.count() << endl; 

		outfile << s << " " << t << endl;
		outfile << Res.size() << endl;
		for(int i = 0; i < Res.size(); i++)
		{
			outfile << Res[i].dis << " " << Res[i].park << endl;
			// cout << Res[i].dis << " " << Res[i].park << endl;
		}
		outfile << endl;

	}
	outfile_res << "GLAD : totol time = "<<GLAD_time_totol.count()<<"\t"<<"average time = "<<GLAD_time_totol.count()/num<<endl;
	outfile_res << Scout << endl;
	outfile_res << "---------------------------------------" << endl;
	outfile_res << endl;

	cout << "GLAD : totol time = "<<GLAD_time_totol.count()<<"\t"<<"average time = "<<GLAD_time_totol.count()/num<<endl;
	cout << Scout << endl;
	

	return 0;
}


