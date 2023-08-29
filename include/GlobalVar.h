
static const int QUEUEMAX = 100000000;

int CELLEDGE =1000;
int NUMOBJECT = 10000;
int CELLSIZE = 1000;
int NUM_THREADS = 40;
int QUERY_SIZE = 1000000;
int K_VALUE = 10;
int ABORTSUM1 = 0;
int ABORTSUM2 = 0;
double RATIO = 0.2;
int DATA = 0; //0: NY, 1:cal, 2:COL
int NUMPOINT = 0;
int ERROR = 0;

struct PA
{
	vector<int> pathnode;
	int dis;
	int park;
	string path;
	PA(vector<int>p,int x,int y){
		pathnode = p;
		dis = x;
		park = y;
		for(int i=0;i<p.size();i++)
		{
		        path+=to_string(p[i]);
		        path+="->";
		}
	}
};

atomic<int> abortNum(0);
atomic<int> queueIndexFront(0);
atomic<int> queueIndexTail(QUERY_SIZE);
vector<int> abortQueryList;
mutex globalMutex;
vector<int> queryList;
struct ObjectInfo {
	int availability = 1;
	mutex locker;
};
vector<ObjectInfo*> OBJECTLIST;

struct INT2{
	int elem[2];
};

