// ChisqTests.cpp -- Chisquare tests on contingency tables
//
// Joe Song
//
// Created: June 16, 2006
// Updated: September 23, 2006.  Added MultinomialChisquareTestTable()
// Modified: September 5, 2008.  Corrected the degrees of freedom
// Last modified: Name changed from multinom.cpp to ChisqTests.cpp
// Modified: September 11, 2011.  Add a new chisq comparison function to make the 
//           comparison more stringent than penalizing by degrees of freedom.
//           Feb 17, 2015. Hua commented all std::out and assert
//#include <iostream>
//using std::cerr;
//using std::cout;
//using std::endl;

#include <cstdlib>
#include <cmath>
#include <cstring>
//#include <cassert>
#include <numeric>
#include <vector>
#include <list>
#include <string>

//Added by Hua Apr 29 2014
#if defined _WIN32
#include <time.h>
#endif
////

using namespace std;

#include <Rcpp.h> // MS Added Feb 8, 2015 to avoid using standard I/O

#include "StatDistributions.h"
#include "ChisqTests.h"
#include "ExactFunctionalTest.h"


bool operator != (const Chisq & c1, const Chisq & c2)
{
    return ! (c1 == c2);
}

bool operator == (const Chisq & c1, const Chisq & c2)
{
    if (c1.m_x != c2.m_x) {
        return false;
    } else if (c1.m_df != c2.m_df) {
        return false;
    } else if (c1.m_pval != c2.m_pval) {
        return false;
    } else {
        return true;
    }
}

vector<int> getRowSums(const vector<vector<int> > & tt)
{
	size_t nRows = tt.size();
	vector<int> rowSums(nRows, 0);
	for(size_t row = 0; row < nRows; row ++) {
		for(size_t col = 0; col < tt[row].size(); col ++) {
			rowSums[row] += tt[row][col];
		}
	}
    
	return rowSums;
}

vector<int> getColSums(const vector<vector<int> > & tt)
{
	size_t nRows = tt.size();
    
    vector<int> colSums;
    
    if(nRows > 0) {
        
        colSums.resize(tt[0].size(), 0);
        
        for(size_t col = 0; col < tt[0].size(); col ++) {
            for(size_t row = 0; row < nRows; row ++) {
                colSums[col] += tt[row][col];
            }
        }
    }
    
	return colSums;
}

int getTotalSum(const vector<vector<int> > & tt)
{
	int totalsum = 0;
	for(size_t row = 0; row < tt.size(); row ++) {
		for(size_t col = 0; col < tt[0].size(); col ++) {
			totalsum += tt[row][col];
		}
	}
    
	return totalsum;
}

double chisq(const vector<vector<int> > & O, const vector<vector<double> > & E)
{
    double chisq = 0.0;
    
    for(size_t i=0; i<O.size(); i++) {
        for(size_t j=0; j<O[i].size(); j++) {
            if(E[i][j] != 0) {
                double d = O[i][j]-E[i][j];
                chisq += d * d / E[i][j];
            }
        }
    }
    
    return chisq;
}

double chisq(const vector<vector<int> > & O)
{
    if (O.size() > 0 ) {
        vector<vector<double> > E = getExpectedTable(O);
        return chisq(O, E);
    } else {
        return 0.0;
    }
}

vector<vector<double> > getExpectedTable(const vector<vector<int> > & contingencytable)
{
    //assert(contingencytable.size()>0);
    vector<vector<double> > expectedtable(contingencytable.size(),
                                          vector<double>(contingencytable[0].size(), 0)
                                          );
    
    vector<int> rowsums = getRowSums(contingencytable);
    vector<int> colsums = getColSums(contingencytable);
    // int totalcounts = getTotalSum(contingencytable);
    int total = accumulate(rowsums.begin(), rowsums.end(), 0);
    for(size_t i=0; i<expectedtable.size(); i++)
    {
        for(size_t j=0; j<expectedtable[i].size(); j++)
        {
            expectedtable[i][j] = (double)rowsums[i]*(double)colsums[j]/(double)total;
        }
    }
    return expectedtable;
}

void ChisquareTest1DNoPValue(const vector<int> & x_obs, 
                             const vector<double> & p_null, 
                             int K, double & chisq, 
                             size_t & df)
{
	int N = 0;

	chisq = 0;
    df = K-1;

	for(int k=0; k<K; k++) {
		N += x_obs[k];
	}
    
	if(N <= 0) {
		return;
	}
    
	for(int k=0; k<K; k++) {
		double x_exp = N * p_null[k];
		if(x_exp != 0) {
			chisq += (x_obs[k] - x_exp)*(x_obs[k] - x_exp)/ x_exp;
		} else if(x_obs[k] != 0) {
			// cerr << "ERROR: expected is zero, but observed is not. Impossible!" << endl;
			// exit(EXIT_FAILURE);
            throw "ERROR: expected is zero, but observed is not. Impossible!";
		}
	}
    
	//legacy code
	//if(N % (int) pow(7.0, 3) == 0) {  // Quick hack to adjust for duplicated sample
	//	chisquare /= pow(7.0, 3);
	//	// cout << ".";
	//}
}

double ChisquareTest(const vector<int> & x_obs, const vector<double> & p_null, int K, double & chisq, size_t & df)
{
    ChisquareTest1DNoPValue(x_obs, p_null, K, chisq, df);
	// return qchisq(K-1, chisq);
    return ChisqPvalue(chisq, K-1);
}


double ChisquareTest(const vector< vector<int> > & table_obs, double & chisq, size_t & df,
                     const vector< vector<double> > & null_prob)
{
	int nrow = (int) table_obs.size();
	int K = (int) table_obs[0].size();
	
    int N = 0;
    
	vector<int> row_sum(nrow), col_sum(K); 

    for(int i = 0; i < nrow; i++) {
    	for(int j = 0; j < K; j++) {
            
			col_sum[j] += table_obs[i][j];
            row_sum[i] += table_obs[i][j]; 
            
		}
		N += row_sum[i];
	}
    	
	chisq = 0;
	
    if(N > 0) {
        for(int i = 0; i < nrow; i++) {
            for(int j = 0; j < K; j++) {
                
                double Eij;
                if (! null_prob.empty()) {
                    Eij = N * null_prob[i][j];
                } else {
                    Eij = row_sum[i] * col_sum[j] / static_cast<double>(N);
                }
                
                if ( Eij > 0 ) {
                    double d = table_obs[i][j] - Eij;
                    
                    chisq += d * d / Eij;
                }                
            }
        }
    }
    
    if( ! null_prob.empty() ) {
        df = nrow * K - 1;  // MS 6/15/2013
    } else {
        df = (nrow-1) * (K-1);
    }
	// return qchisq((int) df, chisq);
    return ChisqPvalue(chisq, (int) df);
}


////////////////
////Hua commented, Feb 15, 2015
////Add by Hua Apr 14 2014, for -M test to do simulation study on contingency table directly.
//#include "TransitionTable.h"
//double applyFunctionalChisqTest(const string & file, int pValMode)
//{
//    TransitionTable tt;
//    
//    tt.scan(file);
//    double p_value;
//    
//	p_value = exact_functional_test(tt.getTransitionTable(), "chisq");
//    
//    /*
//    cout << "pd\tchisqd\tdf" << endl;
//    cout << p_value << "\t" << tt.getChisq()
//    << "\t" << tt.getDf() << endl;
//     */
//    
//    Rcpp::Rcout << "pd\tchisqd\tdf" << endl
//        << p_value << "\t" << tt.getChisq()
//        << "\t" << tt.getDf() << endl;
// 
//	return p_value;
//}
///////////

