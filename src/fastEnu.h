// The fast enumeration algorithm code to obtain P-value
// Enumerate table fit the given row/col sum

#ifndef FASTEUN_H
#define FASTEUN_H

//  fastEnu.h
//    Exact functional test implementaiton using dynamic and
//    quadratic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2023-09-09 (Yiyi Li): Re-write programming to fix bugs
//
//  2021-01-06 (Yiyi Li):Add true hash table
//
//  2019-02-25 (Hien Nguyen): created the namespace DQP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.h of 2.4.7 to EFT_DQP.h to distinguish from other
//     implementations of the exact functional test.

#include "fastEnuNode.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iterator>
#include <unordered_map>

namespace fastEnu
{
    // createKey: convert Rs to equivalent int
    unsigned long int createKey(vector<int> Rs, int layer, int maxCSum);

    // compute the weight between two nodes
    double colChisq(vector<int> & Rs1, vector<int> & Rs2, int sum, const vector<int> &squares, const double &ROWMARGIN);

    double colChisq(vector<int> Rs1, const int &sum, const vector<int> &squares, const double &ROWMARGIN);

    // compute the length from the current node to the end node
    double length(const vector<int> &Rs1, const int &sum, int &layer, const vector<int> &Cs, const vector<double> &factorials);

    // compute the length between two nodes
    double length(const vector<int> &Rs1, const vector<int> &Rs2, const vector<double> &factorials);

    // compute the funchisq without the fixed marginals
    double funchisqByRow(const vector<vector<int>> &observedTable, vector<int> &RSUM, const vector<int> &squares, double &ROWMARGIN);


    // enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
    void createNode(fastEnuNode & node, vector<int> Cs, const vector<int> &Rs, int layer, vector<int> &currRs, int &ncols, int sum1, int sum2,
                    const vector<int> &S, const int &i, const vector<int> &squares, const vector<double> &factorials, vector<fastEnuNode> &Layer,
                    const double &ROWMARGIN, unordered_map<unsigned long int, int> &hashTable, int maxCSum);

    // exact: if true, would only return 0 so that the bound would be evaluate only by DP
    double lower_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN);

    // exact: if ture, sould only return the maximum value. bound would be evaluate only by DP
    double upper_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN);

    // main program for EFT
    double EFTNetwork(const vector<vector<int>> &inputM);
}

#endif
