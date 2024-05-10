//  fastEnu.cpp
//    Exact functional test implementaiton using dynamic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2023-09-09 (Yiyi Li): Re-write programming to fix bugs
//
//  2021-01-06 (Yiyi Li):Add true hash table
//
//  2019-09-07 (Hien Nguyen): fix the error case when the branch and bound occurs at the root node.
//
//  2019-02-25 (Hien Nguyen): insulate code into the namespace DP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.cpp of 2.4.7 to EFT_DQP.cpp to distinguish from other
//     implementations of the exact functional test.
//

#include "fastEnu.h"

using namespace std;

// createKey: convert Rs to equivalent string
unsigned long int fastEnu::createKey(vector<int> Rs, int layer, int maxCSum)
{
  maxCSum++;
  Rs.push_back(layer);
  sort(Rs.begin(), Rs.end());
  unsigned long int eInt = 0;
  for (size_t x = 0; x < Rs.size(); x++)
  {
    eInt *= (unsigned long int)maxCSum;
    eInt += (unsigned long int)Rs[x];
  }
  return eInt;
}

// compute the weight between two nodes
// colchisq = colChisq(network[1][x].getCsum(), sumCol[0], squares, ROWMARGIN);
double fastEnu::colChisq(vector<int> &Rs1, vector<int> &Rs2, int sum, const vector<int> &squares, const double &ROWMARGIN)
{
  if (sum > 0)
  {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs2.size(); x++)
    {
      colchisq += squares[Rs1[x] - Rs2[x]];
    }
    colchisq = colchisq * ROWMARGIN / sum;
    return (colchisq);
  }
  else
    return 0;
}

double fastEnu::colChisq(vector<int> Rs1, const int &sum, const vector<int> &squares, const double &ROWMARGIN)
{
  if (sum > 0)
  {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs1.size(); x++)
    {
      colchisq += squares[Rs1[x]];
    }
    colchisq = colchisq * ROWMARGIN / sum;
    return (colchisq);
  }
  else
    return 0;
}

// compute the length from the current node to the end node
double fastEnu::length(const vector<int> &Rs1, const int &sum, int &layer, const vector<int> &Cs, const vector<double> &factorials)
{
  double length = factorials[sum];
  for (size_t x = 0; x < Rs1.size(); x++)
  {
    length /= factorials[Rs1[x]];
  }
  for (int x = 0; x < layer; x++)
  {
    length /= factorials[Cs[x]];
  }
  return (length);
}

// compute the length between two nodes
double fastEnu::length(const vector<int> &Rs1, const vector<int> &Rs2, const vector<double> &factorials)
{
  double length = 1.0;
  for (size_t x = 0; x < Rs2.size(); x++)
  {
    length /= factorials[Rs1[x] - Rs2[x]];
  }
  return (length);
}

// compute the funchisq without the fixed marginals
double fastEnu::funchisqByRow(const vector<vector<int>> &inputM,
                              vector<int> &RSUM, const vector<int> &squares, double &ROWMARGIN)
{
  double rowchisq = 0.0;
  double funchisq = 0.0;
  for (size_t i = 0; i < inputM.size(); i++)
  {
    rowchisq = 0;
    if (RSUM[i] > 0)
    {
      for (size_t j = 0; j < inputM[0].size(); j++)
      {
        rowchisq += squares[inputM[i][j]];
      }
      rowchisq = rowchisq * ROWMARGIN / RSUM[i];
    }
    funchisq += rowchisq;
  }
  return (funchisq);
}

double fastEnu::lower_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN)
{
  double lower_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q)
  {
    order[q] = q;
  }

  // sort U in increasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2)
       { return U[i1] < U[i2]; });

  for (int l = 0; l < layer; l++)
  {
    // find lower bound for row l
    int runsum = 0;

    for (size_t k = 0; k < nrows; ++k)
    {
      // accumulate the lower bound
      double xavg = (O_colsums[l] - runsum) / (double)(nrows - k);
      if (U[order[k]] < xavg)
      {
        if (O_colsums[l] > 0)
          lower_bound += U[order[k]] * U[order[k]] * ROWMARGIN / (double)O_colsums[l];
        runsum += U[order[k]];
      }
      else
      {
        if (O_colsums[l] > 0)
          lower_bound += (nrows - k) * xavg * xavg * ROWMARGIN / (double)O_colsums[l];
        break;
      }
    }
  }

  return lower_bound;
}

double fastEnu::upper_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN)
{
  double upper_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q)
  {
    order[q] = q;
  }

  // sort U in decreasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2)
       { return U[i1] > U[i2]; });

  for (int l = layer - 1; l >= 0; l--)
  {
    // find lower bound for row l
    if (O_colsums[l] > 0)
    {
      int runsum = 0;

      for (size_t k = 0; k < nrows; ++k)
      {
        // accumulate the lower bound
        int xmax = O_colsums[l] - runsum;
        if (U[order[k]] < xmax)
        {
          if (O_colsums[l] > 0)
            upper_bound += U[order[k]] * U[order[k]] * ROWMARGIN / O_colsums[l];
          runsum += U[order[k]];
        }
        else if (xmax != 0)
        {
          if (O_colsums[l] > 0)
            upper_bound += xmax * xmax * ROWMARGIN / O_colsums[l];
          runsum += xmax;
        }
        else
        {
          break;
        }
      }
    }
  }

  return upper_bound;
}

// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void fastEnu::createNode(fastEnuNode &node, vector<int> Cs, const vector<int> &Rs, int layer, vector<int> &currCs, int &ncols, int sum1, int sum2,
                         const vector<int> &S, const int &i, const vector<int> &squares, const vector<double> &factorials, vector<fastEnuNode> &Layer,
                         const double &ROWMARGIN, unordered_map<unsigned long int, int> &hashTable, int maxCSum)
{
  // When the current column sum completed.
  if (i == ncols)
  {
    // Obtain the length, weight of node
    double len = fastEnu::length(Cs, currCs, factorials);
    int colchisq = fastEnu::colChisq(Cs, currCs, Rs[layer], squares, ROWMARGIN);

    unsigned long int eKey = fastEnu::createKey(currCs, layer, maxCSum);
    unordered_map<unsigned long int, int>::const_iterator keyIt;
    keyIt = hashTable.find(eKey);

    // if the child node does not exist yet, insert it to the next layer as a new node
    if (keyIt == hashTable.end())
    {
      Layer.push_back(fastEnuNode(currCs, eKey));
      node.addChildLink((int)Layer.size() - 1, len, colchisq);
      // update hashTable by insertion
      hashTable.insert(make_pair(eKey, Layer.size() - 1));

      Layer[Layer.size() - 1].setMinPastChisq(node.getMinPastChisq() + colchisq);
      Layer[Layer.size() - 1].setMaxPastChisq(node.getMaxPastChisq() + colchisq);

      Layer[Layer.size() - 1].setLB(fastEnu::lower_bound(layer, currCs, Rs, ROWMARGIN));
      Layer[Layer.size() - 1].setUB(fastEnu::upper_bound(layer, currCs, Rs, ROWMARGIN));

      Layer[Layer.size() - 1].setLengthToEnd(length(currCs, S[layer - 1], layer, Rs, factorials));
    }
    // Changed to fix 0 rowSum problem
    else if (Layer.size() == 0 && keyIt != hashTable.end())
    {
    }
    else
    {
      // if the child node already exists, add a new link from the current node to that child node
      int index = keyIt->second;
      node.addChildLink(index, len, colchisq);

      Layer[index].setMinPastChisq(std::min(Layer[index].getMinPastChisq(), node.getMinPastChisq() + colchisq));
      Layer[index].setMaxPastChisq(std::max(Layer[index].getMaxPastChisq(), node.getMaxPastChisq() + colchisq));
    }
  } // end i == ncols
  else
  {
    // When current column sum not completed, get the next column sum
    int lowerbound, upperbound;

    sum1 += (i > 0 ? Cs[i - 1] : 0);
    sum2 += (i > 0 ? currCs[i - 1] : 0);
    // Get the bound for the value of node.
    lowerbound = std::max(0, Cs[i] - Rs[layer] + sum1 - sum2);
    upperbound = std::min(Cs[i], (layer - 1 >= 0 ? S[layer - 1] : 0) - sum2);

    for (int x = lowerbound; x <= upperbound; x++)
    {
      currCs[i] = x; // Get currecnt column sum at this layer.
      fastEnu::createNode(node, Cs, Rs, layer, currCs, ncols,
                          sum1, sum2, S, i + 1, squares, factorials, Layer, ROWMARGIN, hashTable, maxCSum);
    }
  }
}

double fastEnu::EFTNetwork(const vector<vector<int>> &inputM)
{
  int nrows = (int)inputM.size();
  int ncols = nrows > 0 ? (int)inputM[0].size() : 0;

  if (nrows == 0 || ncols == 0)
  {
    return 1.0;
  }

  int N = 0;
  vector<int> RowSums(nrows, 0);
  vector<int> ColSums(ncols, 0);

  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < ncols; j++)
    {
      N += inputM[i][j];
      RowSums[i] += inputM[i][j];
      ColSums[j] += inputM[i][j];
    }
  }

  if (N == 0)
  {
    return 1;
  }

  int maxCSum = nrows;
  for (int j = 0; j < ncols; j++)
  {
    if (maxCSum < ColSums[j])
    {
      maxCSum = ColSums[j];
    }
  }

  vector<int> squares(N + 1);
  for (int x = 0; x < N + 1; x++)
    squares[x] = x * x;

  vector<double> factorials(N + 1);
  factorials[0] = 1.0;
  for (int x = 1; x <= N; x++)
    factorials[x] = x * factorials[x - 1];
  // for (int x = 1; x <= N; x++)  factorials[x] = factorial<double>(x);

  double marginal = factorials[N];
  for (int x = 0; x < nrows; x++)
  {
    marginal /= factorials[RowSums[x]];
  }
  for (int x = 0; x < ncols; x++)
  {
    marginal /= factorials[ColSums[x]];
  }

  std::vector<int> S(nrows);
  S[0] = RowSums[0];
  for (int x = 1; x < nrows; x++)
  {
    S[x] = S[x - 1] + RowSums[x];
  }

  double ROWMARGIN = 1;
  for (int i = 0; i < nrows; i++)
  {
    if (RowSums[i] > 0)
      ROWMARGIN *= RowSums[i];
  }

  // compute the adjusted funchisq
  double funchisq;
  funchisq = fastEnu::funchisqByRow(inputM, RowSums, squares, ROWMARGIN);

  double funchisqRight = 0;
  for (unsigned i = 0; i < ColSums.size(); i++)
  {
    funchisqRight += squares[ColSums[i]] * ROWMARGIN / N;
  }

  if (funchisq - funchisqRight == 0)
  {
    return 1;
  }

  vector<vector<fastEnuNode>> network(nrows + 1);

  // create the sink node
  vector<int> currCs(ncols, 0);

  network[nrows].push_back(fastEnuNode(ColSums, 0));

  network[nrows][0].addPastLen(1.0, 0);
  network[nrows][0].setMaxPastChisq(0);
  network[nrows][0].setMinPastChisq(0);

  network[nrows][0].setLB(fastEnu::lower_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setUB(fastEnu::upper_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setLengthToEnd(fastEnu::length(ColSums, S[nrows - 1], nrows, RowSums, factorials));

  if (network[nrows][0].getLB() >= funchisq)
  {
    return 1;
  }

  // hash table
  unordered_map<unsigned long int, int> hashTable;

  // generate the network

  for (int layer = nrows; layer > 1; layer--)
  {
    for (size_t n = 0; n < network[layer].size(); n++)
    {
      if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq)
      {
        if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq)
        {
          fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
        }
      }
    }
  }

  double rowchisq;
  for (size_t x = 0; x < network[1].size(); x++)
  {
    rowchisq = colChisq(network[1][x].getCsum(), RowSums[0], squares, ROWMARGIN);
    network[1][x].setLB(rowchisq);
    network[1][x].setUB(rowchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= nrows; layer++)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {

      // network[layer][node].setLengthToEnd(length(network[layer][node].getCsum(), S[layer - 1], layer, RowSums, factorials));

      minLB = DBL_MAX;
      maxUB = 0;

      if (network[layer][node].getSize() > 0)
      {

        for (int child = 0; child < network[layer][node].getSize(); child++)
        {

          rowchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = rowchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound)
            minLB = tempBound;

          tempBound = rowchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound)
            maxUB = tempBound;
        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }
    }
  }

  double lengthSoFar = 1;
  double chisqSoFar = 0;
  double pvalue = 0;
  int pastSize = 0;

  // traverse the network in a breadth-first strategy
  for (int layer = nrows; layer >= 1; layer--)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++)
      {

        chisqSoFar = network[layer][node].getPastChisq(i);

        if (
            chisqSoFar + network[layer][node].getUB() < funchisq)
        {
        }
        else
        {

          lengthSoFar = network[layer][node].getPastLen(i);

          if (chisqSoFar + network[layer][node].getLB() >= funchisq)
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else
          {
            for (int k = 0; k < network[layer][node].getSize(); k++)
            {
              network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                  network[layer][node].getLengthToChildren(k) * lengthSoFar,
                  network[layer][node].getColChisqToChildren(k) + chisqSoFar);
            } // end child
          }
        }
      } // end pastChisq

    } // end node
  }   // end layer
  hashTable.clear();
  pvalue /= marginal;

  return pvalue;
}
