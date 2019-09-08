
#ifndef NODE_H
#define NODE_H

#include <vector>

using namespace std;

class Node {

public:
	Node();
	Node(vector<int> Rs);
	Node(vector<int> Rs, int eInt);

	void show();

	void addRsum (int Rs);
	vector<int> getRsum();

	int getChildrenIndex(int x);
	int isChildInList(int x);

	void addLength(double x, int index);
	double getLengthToChildren(int x);

	void setColChisqToChildren(int x,  int colchisq);
	int getColChisqToChildren(int x);

	void addPastLen(double x,  int y);
	double getPastLen(int x);

	//void addPastChisq(int x);
	int getPastChisq(int x);
	int bSearch(int chisq);
	int getPastChisqSize();

	int getPastSize();

	void setUB(unsigned int x);
	unsigned int getUB();

	void setLB(unsigned int x);
	unsigned int getLB();

	void setEquiInt(int x);
	int getEquiInt();

	int getSize();

	void setLengthToEnd(double x);
	double getLengthToEnd();

	void quicksort(int left, int right);

	void addChildLink(int index, double len, int colchisq);

	void setMinPastChisq( int x);
	 int getMinPastChisq();

	void setMaxPastChisq( int x);
	 int getMaxPastChisq();

private:
	vector<int> Rsum;	// the remaining rowsums after the previous columns are enumerated
	int equiInt;		// the equivalent integer converted from Rsum to compare the nodes faster
	double lengthToEnd; // the length from this node to the end node, used when the whole branch is counted

	unsigned int ub;			// the upper bound for this node
	unsigned int lb;			// the lower bound for this node

	vector<int> ChildrenIndex;			// to store the indices of the children node
	vector<double> lengthToChildren;	// to store the lenghths to the children node
	vector< int> colChisqToChildren;	// to store the weights to the children node (weight = partial funchisq for the column enumerated)

	vector<double> pastLen;		// the cummulative lengths from the start node to this node, each entry of this list will be called lengthSoFar
	vector<int> pastChisq;	// the cummulative weights from the start node to this node, each entry of this list will be called chisqSoFar or weightSoFar

	vector<vector<pair <int, int>>> nodeTable;

	int minPastChisq;
	int maxPastChisq;
};

#endif
