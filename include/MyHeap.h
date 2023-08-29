/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/

#ifndef MYHEAP_H
#define MYHEAP_H
#include<vector>

using namespace std;
class MyHeap
{
public:
	struct heapElement {
		int element;
		int priority;
	};
	vector<int> position;
	vector<heapElement> heap;
	void clear();
	void init(int _numElement);
	void push(int _element, int _priority);
	void update(int _element, int _newPriority);
	void pop();
	MyHeap::heapElement top();
	bool empty();
	bool hasElement(int _element);
	int getPriority(int _element);
	void rise(int _index);
	void sink(int _index);
	int parent(int _index);
	int child(int _index);
	void swap(int _index1, int _index2);
	MyHeap();
	~MyHeap();
};
#endif
