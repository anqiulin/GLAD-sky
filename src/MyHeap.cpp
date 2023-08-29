/*
	Author: Dan He
	Email: d.he@uq.edu.au
	Date: 14/08/2018
*/

#include "MyHeap.h"

void MyHeap::init(int _numElement) {
	position.assign(_numElement, -1);
}

void MyHeap::push(int _element, int _priority) {
	int eleIndex = heap.size();
	position[_element] = eleIndex;
	heapElement item;
	item.element = _element;
	item.priority = _priority;
	heap.push_back(item);
	rise(eleIndex);
}

void MyHeap::update(int _element, int _newPriority) {
	int eleIndex = position[_element];
	heap[eleIndex].priority = _newPriority;
	rise(eleIndex);
}

void MyHeap::pop() {
	position[heap[0].element] = -1;
	if (heap.size() > 1) {
		heap[0] = heap[heap.size() - 1];
		position[heap[0].element] = 0;
	}
	heap.pop_back();
	if (!heap.empty()) {
		sink(0);
	}
}

MyHeap::heapElement MyHeap::top() {
	return heap[0];
}

bool MyHeap::empty() {
	return heap.empty();
}

bool MyHeap::hasElement(int _element)
{
	return position[_element] != -1;
}

int MyHeap::getPriority(int _element) {
	int index = position[_element];
	return heap[index].priority;
}

void MyHeap::rise(int _index) {
	while (_index > 0 && heap[parent(_index)].priority > heap[_index].priority)
	{
		swap(_index, parent(_index));
		_index = parent(_index);
	}
}

void MyHeap::sink(int _index) {
	while (true) {
		int child1 = child(_index), child2 = child1 + 1;
		int min = _index;
		if (child1 < heap.size() && heap[child1].priority < heap[min].priority)
			min = child1;
		if (child2 < heap.size() && heap[child2].priority < heap[min].priority)
			min = child2;
		if (min != _index)
		{
			swap(_index, min);
			_index = min;
		}
		else
			return;
	}
}

int MyHeap::parent(int _index) {
	return (_index - 1) / 2;
}

int MyHeap::child(int _index) {
	return _index*2 + 1;
}

void MyHeap::swap(int _index1, int _index2) {
	position[heap[_index1].element] = _index2;
	position[heap[_index2].element] = _index1;
	std::swap(heap[_index1], heap[_index2]);
}

MyHeap::MyHeap()
{
}

void MyHeap::clear() {
	for (auto &item : heap)
		position[item.element] = -1;
	heap.clear();
}

MyHeap::~MyHeap()
{
}
