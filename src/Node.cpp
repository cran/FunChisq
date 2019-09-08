#include "Node.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

Node::Node()
{

}

Node::Node(vector<int> Rs, int eInt) {
	Rsum = Rs;
	equiInt = eInt;
	nodeTable.resize(199);
}

Node::Node(vector<int> Rs)
{
	Rsum = Rs;
	std::sort(Rs.begin(), Rs.end());

	int eInt = 0;
	for (size_t x = 0; x < Rs.size(); x++) {
		eInt *= 127;
		eInt += Rs[x];
	}
	equiInt = eInt;

}

void Node::show()
{
  /*
	for (size_t i = 0; i < Rsum.size(); i++) {
		cout << Rsum[i] << " ";
	}
	cout << " :: ";
	for (size_t i = 0; i < ChildrenIndex.size(); i++) {
		cout << ChildrenIndex[i] << " ";
	}
	cout << "::" << ub << " :: " << lb << endl;
  */
}

void Node::addRsum(int Rs)
{
	Rsum.push_back(Rs);
}

vector<int> Node::getRsum()
{
	return Rsum;
}

int Node::getChildrenIndex(int x)
{
	return ChildrenIndex[x];
}

double Node::getLengthToChildren(int x) {
	return lengthToChildren[x];
}

int Node::getColChisqToChildren(int x) {
	return colChisqToChildren[x];
}

void Node::setUB(unsigned int x) {
	ub = x;
}
void Node::setLB(unsigned int x) {
	lb = x;
}

unsigned int Node::getLB() {
	return lb;
}

unsigned int Node::getUB() {
	return ub;
}

void Node::setEquiInt(int x) {
	equiInt = x;
}

int Node::getEquiInt() {
	return equiInt;
}

int Node::getSize() {
	return ChildrenIndex.size();
}

void Node::addPastLen(double x, int y) {


	int i = y % nodeTable.size();
	size_t j = 0;
	while (j < nodeTable[i].size() && nodeTable[i][j].first != y) j++;

	if (j < nodeTable[i].size()) pastLen[nodeTable[i][j].second] += x;

	else {
		pastLen.push_back(x);
		pastChisq.push_back(y);

		//update hash table:
		nodeTable[i].push_back(std::make_pair(y, pastChisq.size()-1));
	}

	/*
	int i = 0;
	while (i < pastChisq.size() && (pastChisq[i] != y)) i++;

	if (i < pastChisq.size()) {
		pastLen[i] += x;
	}
	else {
		pastLen.push_back(x);
		pastChisq.push_back(y);
	}
	*/
}

double Node::getPastLen(int x) {
	return pastLen[x];
}

int Node::getPastSize() {
	return pastLen.size();
}

int Node::getPastChisq(int x) {
	return pastChisq[x];
}

void Node::setLengthToEnd(double x) {
	lengthToEnd = x;
}

double Node::getLengthToEnd() {
	return lengthToEnd;
}

int Node::getPastChisqSize() {
	return pastChisq.size();
}

void Node::quicksort(int left, int right) {

	double pivot = pastChisq[(left + right) / 2];
	int i, j;
	//double temp;
	i = left;
	j = right;

	while (i <= j) {
		while (pastChisq[i] < pivot) i++;
		while (pastChisq[j] > pivot) j--;

		if (i <= j) {
			//swap arr[i] and arr[j]
			std::swap(pastChisq[i], pastChisq[j]);
			std::swap(pastLen[i], pastLen[j]);
			i++;
			j--;
		}
	}

	if (left < j) quicksort(left, j);
	if (i < right) quicksort(i, right);
}


int Node::bSearch(int chisq) {
	auto it = std::lower_bound(pastChisq.begin(), pastChisq.end(), chisq);
	std::size_t index = std::distance(pastChisq.begin(), it);
	return index;
}


int Node::isChildInList(int x) {
	size_t i = 0;
	while (i < ChildrenIndex.size() && x != ChildrenIndex[i]) i++;
	if (i < ChildrenIndex.size()) 	return i;
	return -1;
}

void Node::addLength(double x, int index) {
	lengthToChildren[index] += x;
}


void Node::addChildLink(int index, double len, int colchisq) {
	ChildrenIndex.push_back(index);
	lengthToChildren.push_back(len);
	colChisqToChildren.push_back(colchisq);
}

void Node::setColChisqToChildren(int x, int colchisq) {
	colChisqToChildren[x] = colchisq;
}


void Node::setMinPastChisq(int x) {
  minPastChisq = x;
}
int Node::getMinPastChisq() {
  return (minPastChisq);
}

void Node::setMaxPastChisq(int x) {
  maxPastChisq = x;
}
int Node::getMaxPastChisq() {
  return maxPastChisq;
}
