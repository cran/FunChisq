//
//  main.cpp
//  eft-mhg
//
//  Created by Hua Zhong on 6/28/15.
//  Copyright (c) 2015 New Mexico State University. All rights reserved.
//

#include "define.h"
#include <vector>
//#include "ExactFunctionalTest.h"//Commented by Hua, Jun 17 2015
#include <Rcpp.h>
//#include <classic/RcppMatrix.h>
#include "ExactFunctionalTest.h"
#include "StatDistributions.h"

//using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double ExactFunctionalTest(const IntegerMatrix & nm) {
  // RcppExport SEXP ExactFunctionalTest(SEXP x){// x is a numeric matrix
  //RcppMatrix<int> nm(x);
  //vector< vector< int > > C = nm.stlMatrix();
  //double pVal = exact_functional_test(C, "chisq");//Commented by Hua, Jun 17 2015

  // The following is added to replace RcppMatirx to remove dependence on RcppClassic:
  // IntegerMatrix nm(x);

  vector< vector< int > > C(nm.nrow(), vector<int>(nm.ncol()));
  for(size_t i = 0; i < nm.nrow(); ++i) {
    for(size_t j = 0; j < nm.ncol(); ++j) {
      C[i][j] = nm(i, j);
    }
  }

  mydouble fc = 0.0;
  mydouble pVal = exact_func_test_multi_hypergeometric
    (C, fc, LBON, UB_BY_ROW, (enum PVAL) PVAL);

  //     mydouble pVal = exact_func_test_multi_hypergeometric
  //       (C, fc, LBOFF, UBOFF, (enum PVAL) PVAL);

  //return wrap(pVal);
  return pVal;
}

// [[Rcpp::export]]
DataFrame interactions(const IntegerMatrix & expression_matrix, const List & parent_index,
                       const IntegerVector & child_index, const String & index_kind) {

  //Take inputs
  vector< vector< int > > x (expression_matrix.nrow(), vector<int> (expression_matrix.ncol(), 0));
  for(size_t i = 0; i < expression_matrix.nrow(); ++i) {
    for(size_t j = 0; j < expression_matrix.ncol(); ++j) {
      x[i][j] = expression_matrix(i, j);
    }
  }

  vector< vector< int > > P;
  for(size_t i = 0; i < parent_index.length(); ++i) {
    P.push_back(as< vector < int > >(parent_index[i]));
  }

  vector< int > C = as< vector < int > >(child_index);

  string method = index_kind;
  ////

  vector < mydouble > p_value (C.size(), 1);
  vector < mydouble > statistic (C.size(), 0);
  vector < mydouble > estimate (C.size(), 0);

  for (size_t i=0; i<C.size(); i++) {
    vector<int> P_index = P[i]; //1st index is 0
    for(size_t j=0; j<P_index.size();j++){
      P_index[j] = P_index[j] - 1;
    }
    int C_index = C[i] - 1;

    vector< vector < int > > P_exp;
    for(size_t j=0; j<P_index.size();j++){
      P_exp.push_back(x[P_index[j]]);
    }
    vector <int> C_exp = x[C_index];

    //minimun level require to be 0;
    //parent max level for each parent
    vector<int> P_max_level(P_index.size(), 0);
    for(size_t j=0;j<P_index.size();j++){
      P_max_level[j] = *max_element(P_exp[j].begin(), P_exp[j].end());
    }
    //chile max level
    int C_max_level = *max_element(C_exp.begin(), C_exp.end());

    //Build contingency table
    //size: (((max_P_1+1) * (max_P_2+1) * ... * (max_P_i+1)) * (C_max_level+1)
    //Parent combinatorial row:
    //    0   0   0
    //    0   0   1
    //    0   1   0
    //    0   1   1
    //    1   0   0
    //    ...
    //    1   1   1
    int T_P_levels = 1;
    for(size_t j=0;j<P_index.size();j++){
      T_P_levels *= (P_max_level[j] + 1);
    }

    vector < vector < int > > T (T_P_levels, vector <int> (C_max_level+1, 0));

    for(size_t j=0; j<C_exp.size(); j++){
      int accumulation_tmp = 0; // table combined parent index
      int product_tmp = 1;
      for(int k=(int)P_index.size()-1; k>=0; k--){
        if(k!=P_index.size()-1) product_tmp *= (P_max_level[k+1]+1);
        accumulation_tmp += P_exp[k][j] * product_tmp;
      }

      T[accumulation_tmp][C_exp[j]]++;
    }

    //Remove all-0 rows and columns.
    //row
    size_t j=0;
    while(j < T.size()){
      int row_all_zero = true;
      for(size_t k=0; k<T[j].size(); k++){
        if(T[j][k]!=0) {
          row_all_zero = false;
          break;
        }
      }
      if(row_all_zero == true){
        T.erase(T.begin() + j);
        continue;
      }
      j++;
    }
    //column
    j=0;
    while(j < T[0].size()){
      int column_all_zero = true;
      for(size_t k=0; k<T.size(); k++){
        if(T[k][j]!=0) {
          column_all_zero = false;
          break;
        }
      }
      if(column_all_zero == true){
        for(size_t k=0; k<T.size(); k++){
          T[k].erase(T[k].begin() + j);
        }
        continue;
      }
      j++;
    }

    //calculating statistics
    mydouble estimate_tmp = 0;
    mydouble statistic_tmp = funchisq(T, estimate_tmp, method);
    mydouble p_value_tmp = ChisqPvalue(statistic_tmp, ((int)T.size()-1) * ((int)T[0].size()-1));

    p_value[i] = p_value_tmp;
    statistic[i] = statistic_tmp;
    estimate[i] = estimate_tmp;
  }

  DataFrame output = DataFrame::create(Named("p.value")=p_value,
                                       Named("statistic")=statistic,
                                       Named("estimate")=estimate);
  return output;
}
