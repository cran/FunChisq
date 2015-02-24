// TransitionTable.cpp
//
// Joe Song
// Created: September 27, 2008
// Modified:
//   September 12, 2011. Added different options for parent comparison
//   January 24, 2013.  Extracted friend functions for comparative chisq
//     analysis to ComparativeChisq.cpp
//  Feb 18, 2015. Hua commented all std::~, and assert

//#include <iostream>
//#include <fstream>
//#include <sstream>
//using std::cerr;
//using std::cout;
//using std::endl;
//using std::string;
//using std::stringstream;
//using std::ifstream;
//using std::ofstream;
////using std::exit;

#include <algorithm>
//#include <cassert>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "TransitionTable.h"
#include "ChisqTests.h"

TransitionTable operator - (const TransitionTable & tt1, const TransitionTable & tt2)
{
	TransitionTable tt = tt1;
    
	for(size_t r=0; r < tt1.m_transitionTable.size(); r++) {
		for(size_t c=0; c < tt1.m_transitionTable[r].size(); c++) {
			tt.m_transitionTable[r][c] = (tt1.m_transitionTable[r][c] - tt2.m_transitionTable[r][c]); // absolute value?
		}
	}
	
	return tt;
}

TransitionTable & abs(TransitionTable & tt)
{
	for(size_t r=0; r < tt.m_transitionTable.size(); r++) {
		for(size_t c=0; c < tt.m_transitionTable[r].size(); c++) {
			tt.m_transitionTable[r][c] = abs(tt.m_transitionTable[r][c]);
		}
	}
	return tt;
}


/*There is a potential off-by-one error when using this function.
It assumes IDs stored in the incoming vector 'parents' start from 1
*/

const TransitionTable &
TransitionTable::operator += (const TransitionTable & tt)
{
	// check compatibility of the two tables
	if( m_transitionTable.size() > 0 ) {

        if(m_transitionTable.size() != tt.m_transitionTable.size() ||
           m_transitionTable[0].size() != tt.m_transitionTable[0].size()) {
            
            /* MS Commented Feb 8, 2015
            cerr << "ERRRO in TransitionTable::operator +=: incompatible tables!"
            << endl;
            
            exit(EXIT_FAILURE);
            */
            throw "ERRRO in TransitionTable::operator +=: incompatible tables!";
        }
        
		for(size_t r = 0; r < nrow(); r++) {
			for(size_t c = 0; c < m_transitionTable[r].size(); c++) {
				m_transitionTable[r][c] += tt.m_transitionTable[r][c];
			}
		}
        
	} else {
        
        *this = tt;
        
    }
    
	return *this;
}


vector<int> TransitionTable::getRowSums() const
{
    return ::getRowSums(m_transitionTable);
}

//added by Yang Zhang 1.2.2013
vector<int> TransitionTable::getColSums() const
{
    return ::getColSums(m_transitionTable);
}

int TransitionTable::getTotalSum() const
{
    return ::getTotalSum(m_transitionTable);
}


void TransitionTable::reset()
{
    for(size_t i = 0; i < m_transitionTable.size(); i++)
        for(size_t j = 0; j < m_transitionTable[i].size(); j++)
            m_transitionTable[i][j] = 0;
}
