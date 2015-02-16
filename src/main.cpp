#include <vector>
#include "ExactFunctionalTest.h"
#include <Rcpp.h>
#include <classic/RcppMatrix.h>

using namespace Rcpp;

RcppExport SEXP ExactFunctionalTest(SEXP x){// x is a numeric matrix
	RcppMatrix<int> nm(x);
	vector< vector< int > > C = nm.stlMatrix();
	double pVal = exact_functional_test(C, "chisq");
	return (wrap(pVal));
}
