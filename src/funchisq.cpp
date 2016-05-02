// funchisq.cpp
//
// Created: April 30, 2016. Extracted from ExactFunctionalTest.cpp

#include "define.h"
#include <vector>
#include <string>
#include <cmath>

mydouble funchisq(const std::vector<std::vector<int> > & O, mydouble & estimate,
                  const std::string index_kind){
  mydouble fc = 0.0;
  if (O.size() == 0) {
    return fc;
  } else if(O[0].size() == 0) {
    return fc;
  }

  int n = 0;
  std::vector <int> colsums ((int)O[0].size(), 0);
  std::vector <int> rowsums ((int)O.size(), 0);

  for (size_t i=0; i<O.size(); i++) {
    for (size_t j=0; j<O[i].size(); j++) {
      n += O[i][j];
      colsums[j] += O[i][j];
      rowsums[i] += O[i][j];
    }
  }

  if(n == 0)return fc;

  size_t nrows = O.size();  // number of rows
  size_t ncols = O[0].size();  // number of columns

  mydouble ej = n / (mydouble) ncols;
  mydouble col_chisq = 0.0;

  if(ej>0){
    for (size_t j=0; j<ncols; ++j) {
      col_chisq += (colsums[j] - ej) * (colsums[j] - ej) / ej;
    }
  }
  fc -= col_chisq;

  for (size_t i=0; i<nrows; ++i) {
    // Expected cound for cell (i,j):
    mydouble eij = rowsums[i] / (mydouble) ncols;
    if (eij > 0) {
      for (size_t j=0; j<ncols; ++j) {
        fc += (O[i][j] - eij) * (O[i][j] - eij) / eij;
      }
    }
  }

  if (index_kind == "conditional") {
    estimate = std::sqrt(std::abs(fc) / (n * (ncols - 1) - col_chisq));
  }else if(index_kind == "unconditional"){
    estimate = std::sqrt(std::abs(fc) / (n * (ncols - 1)));
  }

  return fc;
}

mydouble funchisq(const std::vector<std::vector<int> > & O, const std::vector<int> & rowsums,
                  const std::vector<int> & colsums, int n)
{
  mydouble fc = 0.0;

  if (n == 0 || O.size() == 0) {
    return fc;
  } else if(O[0].size() == 0) {
    return fc;
  }

  size_t nrows = O.size();  // number of rows
  size_t ncols = O[0].size();  // number of columns

  mydouble ej = n / (mydouble) ncols;
  if(ej>0){
    for (size_t j=0; j<ncols; ++j) {
      fc -= (colsums[j] - ej) * (colsums[j] - ej) / ej;
    }
  }

  for (size_t i=0; i<nrows; ++i) {
    // Expected cound for cell (i,j):
    mydouble eij = rowsums[i] / (mydouble) ncols;
    if (eij > 0) {
      for (size_t j=0; j<ncols; ++j) {
        fc += (O[i][j] - eij) * (O[i][j] - eij) / eij;
      }
    }
  }
  return fc;
}
