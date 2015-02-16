// TransitionTable.h
//
// Joe Song
// Created: September 27, 2008

#pragma once

#include <vector>
using std::vector;

#include <iostream>
#include <fstream>
using namespace std;

// History:
//   Added by Haizhou Wang, Feb 13, 2013
//   #define INITIAL_P_VALUE 1.0 // MS 2/23/2013 Changed from 1.01 to 1.0
//   MS 6/1/2013 Changed to the following 
const double INITIAL_P_VALUE = 1.01;

class TransitionTable
{
	// friend bool compare(const TransitionTable & tt1, const TransitionTable & tt2);
        
	friend TransitionTable operator - (const TransitionTable & tt1, const TransitionTable & tt2);
	friend TransitionTable & abs(TransitionTable & tt);
    
public:
    
	vector<int> getRowSums() const;

	//added by Yang Zhang 2013.1.2
	vector<int> getColSums() const;
	int getTotalSum() const;
    
	const TransitionTable & operator += (const TransitionTable & tt);
    
	int get(int r, int c) const { return m_transitionTable[r][c]; }
	void set(int r, int c, int val) { m_transitionTable[r][c] = val; }
    
	const vector< vector<int> > & getTransitionTable() const { return m_transitionTable; }
    
	void setTransitionTable(vector< vector<int> > transtable) { m_transitionTable = transtable; }
    
	void setpValue(const double pvalue) { m_pValue = pvalue; }
	void setAdjustedpValue(const double pvalue) { m_adjustedpValue = pvalue; }
    
	void setChisq(const double chisq) { m_chisq = chisq; }
    
	void setDf(const int df) { m_df = df; }
    
	//added by yangzhang 11.26.2008
    
	double getpValue() const { return m_pValue; }
	double getAdjustedpValue() const { return m_adjustedpValue; }
    
	double getChisq() const { return m_chisq; }
    
	int getDf() const { return m_df; }

    // added by Tyler Hunt 08/01.2012 //
    // reset the table to contain only zeros
    void reset();
    
    bool scan(const string & file);

    size_t nrow () const { return m_transitionTable.size(); }
    size_t ncol () const {
        if(m_transitionTable.size()>0) {
            return m_transitionTable[0].size();
        } else {
            return 0;
        }
    }

    
protected:
    /*
     Comment added by Haizhou Wang, July 23, 2012
     This m_transitionTable actually is a contingency table.
     It records how many times each child value show up under different parent combination
     */
	vector< vector<int> > m_transitionTable;
    
	double m_chisq;		// chisquare value
	int m_df;			// degrees of freedom
	double m_pValue;	// p-value
    
	double m_adjustedpValue; // p-value adjusted for multiple comparison
    
	//int m_pValueMode; // chisquare test method
    
	// vector<int> m_truthtable;
    
};

