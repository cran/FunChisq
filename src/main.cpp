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
