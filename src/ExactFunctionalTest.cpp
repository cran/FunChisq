//
//  EMTFisherTests.cpp
//  gln
//
//  Created by Joe Song on 2/12/13.
//
//
#include <cmath>
#include <ctime>
#include <numeric>

using namespace std;

#include "ExactFunctionalTest.h"

double factorial_vector(vector<int>);//Hua Added, May 2 2014

//----------------------------------------------------------------------------

//compute fisher's probability
double FisherProb(const vector<vector<int> > & A)
{
    // Input:
    //   A -- contingency table
    
	// double product = 1.0; Commented by MS Feb 8, 2015
    
    //Modified by Hua, May 2 2014
	vector<int> rowsums = getRowSums(A);
	vector<int>	colsums = getColSums(A);
	vector<int> v;
	v.reserve(rowsums.size() + colsums.size());
	v.insert(v.end(), rowsums.begin(), rowsums.end());
	v.insert(v.end(), colsums.begin(), colsums.end());
	int total = accumulate(rowsums.begin(), rowsums.end(), 0);
	v.push_back(-total);
    
	for(size_t i = 0; i<A.size(); i++)
	{
		for(size_t j=0; j<A[i].size(); j++)
		{
			//product *= factorial(A[i][j]);//Hua comment it, May 2 2014
			v.push_back(-A[i][j]);//Added by Hua, May 2 2014
		}
	}
    
	return factorial_vector(v);//Added by Hua, May 2 2014
    ////

}

void EMTFunctionalChisq::initialize_customized_row_col_sum(const vector<TransitionTable> & Cs){
		m_observedChisq.resize(Cs.size());
        m_nullChisq.resize(Cs.size());
        m_colSumChisq.resize(Cs.size());
        m_boundChisqs.resize(Cs.size());
        //        m_boundHeteroChisq = 0;
        //        m_row=0;
        
		vector<double> moreOrLessExtreme_tmp;//Hua added, May 17 2014
		moreOrLessExtreme_tmp.resize(Cs.size());//Hua added, May 17 2014
        
        for (size_t k=0; k<Cs.size(); k++) {
            double chisq;
            size_t df;
            double p_tmp = ChisqDirTest(Cs[k].getTransitionTable(), chisq, df);
            m_observedChisq[k] = chisq;
            
			//Hua added, May 17 2014
            moreOrLessExtreme_tmp[k] = p_tmp;
			////
            
            m_totalObservedChisq += m_observedChisq[k];
        }
        
        //Hua added, Jun 5 2014
        for (size_t k=0; k<Cs.size(); k++) {
            double min=moreOrLessExtreme_tmp[0];
            if(min < 1-moreOrLessExtreme_tmp[0])min = 1-moreOrLessExtreme_tmp[0];
            int minIndex=0;
            int index=0;
            while(index < moreOrLessExtreme_tmp.size()){
                if(moreOrLessExtreme_tmp[index] < min || 1-moreOrLessExtreme_tmp[index] < min){
                    minIndex = index;
                }
                index++;
            }
            if(moreOrLessExtreme_tmp[minIndex]<=0.5){
                m_moreOrLessExtreme = false;
            }else{
                m_moreOrLessExtreme = true;
            }
        }
        ////
        
        ////Add by Hua Mar 20 2014
        for (size_t k=0; k<Cs.size(); k++) {
            size_t df=0;
            int nrow = (int) Cs[k].getTransitionTable().size();
            int K = (int) Cs[k].getTransitionTable()[0].size();
            
            vector<double> p_null(K, 1.0/K);
            vector<int> n_y(K);  // the histogram of child y
            
            for(int j = 0; j < K; j++){
                
                for(int i = 0; i < nrow; i++){
                    
                    n_y[j] += Cs[k].getTransitionTable()[i][j];
                    
                }//end for
                
            }//end for
            
            double chisq_y;
            ChisquareTest1DNoPValue(n_y, p_null, K, chisq_y, df);
            m_colSumChisq[k] = chisq_y;
        }
        /////
}

void EMTFunctionalChisq::initialize(const vector<TransitionTable> & Cs)
    {
        // compute required row sums
		m_requiredRowSums.resize(Cs.size());
		for(size_t k=0; k < Cs.size(); k++) {
			m_requiredRowSums[k] = Cs[k].getRowSums();
		}
    
		// compute required column sums
		m_requiredColSums.resize(Cs.size());
		for(size_t k=0; k < Cs.size(); k++) {
			m_requiredColSums[k] = Cs[k].getColSums();
		}

        initialize_customized_row_col_sum(Cs);
    }
    
void EMTFunctionalChisq::processTable(size_t k, const EMTEnumerator & e,
                              const vector<TransitionTable> & Cs)
    {
        double chisq;
        size_t df;
        ChisqDirTest(e.As[k].getTransitionTable(), chisq, df);
        m_nullChisq[k] =chisq;
    }
    
bool EMTFunctionalChisq::isMoreExtreme() const
    {
        double nullTotalChisq = 0, observedTotalChisq = 0;
        
        for(size_t k=0; k < m_nullChisq.size(); k++) {
            nullTotalChisq += m_nullChisq[k];
            observedTotalChisq += m_observedChisq[k];
        }
        
		if(m_moreOrLessExtreme ==false){//Hua added, May 17 2014
			return nullTotalChisq >= observedTotalChisq;
		}else{
            return nullTotalChisq < observedTotalChisq;//No equal, modified by Hua, Jun 5 2014
		}
    }

double EMTFunctionalChisq::evaluate(const EMTEnumerator & e, const vector<TransitionTable> & Cs)
    {
        // Multiple independent Fisher's exact tests
        double P;
        
		vector<vector<int> > row = e.ARowsums;
		vector<vector<int> > col = e.AColsums;
		
        if( isMoreExtreme() ) {
            P = 1;
            for(size_t k=0; k < Cs.size(); k++) {
				P *= FisherProb(e.As[k].getTransitionTable());
            }
        } else {
            P = 0;
        }
        
        return P;
    }
    
	vector<TransitionTable> EMTFunctionalChisq::generateTables(const vector<TransitionTable> & Cs) const
	{
		vector<TransitionTable> As(Cs);
    
		for (size_t k=0; k<Cs.size(); k++) {
			As[k].reset();
		}
    
		return As;
	}


    //Added by Hua Mar 20 2014
	string EMTFunctionalChisq::bound(size_t k, size_t i, size_t j, const EMTEnumerator & e, const vector<TransitionTable> & Cs){
		string result = "not-to-skip"; // "to keep entire branch"
        
        int rowNum = e.As[k].getTransitionTable().size();
        int colNum = e.As[k].getTransitionTable()[0].size();
        
		
        
        if(i>0 && i<rowNum && j==0){
			m_skip = false;
			double boundHeteroChisq = 0;//m_boundHeteroChisq;
        
			for (size_t r=0; r<i; r++) {
        		double chisq_tmp;
				size_t df=0;
				vector<double> p_null(colNum, 1.0/colNum);
				ChisquareTest1DNoPValue(e.As[k].getTransitionTable()[r], p_null, colNum, chisq_tmp, df);
		 		boundHeteroChisq += chisq_tmp;
			}
        
			vector<int> colTotalSum = m_requiredColSums[k];//Cs[k].getColSums();
        
			for (size_t c=0; c<colTotalSum.size(); c++) {
        		colTotalSum[c] = colTotalSum[c] - e.AColsums[k][c];
			}
        
			sort(colTotalSum.begin(), colTotalSum.end());
        
			vector<int> rowTotalSum = m_requiredRowSums[k];//s[k].getRowSums();
        
			vector<int> tmpRow;
        
			if(m_moreOrLessExtreme == false){//Hua added, May 17 2014
        		for (size_t r=i; r <rowNum; r++) {
        			tmpRow.resize(colTotalSum.size());
        			int balls = rowTotalSum[r];
        			int index = colNum-1;
        			while (balls > 0) {
        				if(balls <= colTotalSum[index]){
        					tmpRow[index] = balls;
        					balls = 0;
        				}else{
        					tmpRow[index] = colTotalSum[index];
        					balls -= colTotalSum[index];
        					index--;
        				}
        			}
        
        			double chisq_tmp;
        			size_t df=0;
        			vector<double> p_null(colNum, 1.0/colNum);
        			ChisquareTest1DNoPValue(tmpRow, p_null, colNum, chisq_tmp, df);
        			boundHeteroChisq += chisq_tmp;
        			tmpRow.clear();
        		}
        
        		result = (boundHeteroChisq - m_colSumChisq[k] +1e-07 < m_totalObservedChisq) ?
        		"to skip entire branch" : "not-to-skip";
			}else{
        		//Hua added, May 17 2014
        		// compare lowest FunChisq and observed FunChisq
        		for (size_t r=i; r <rowNum; r++) {
        			tmpRow.resize(colTotalSum.size());
        			int balls = rowTotalSum[r];
        			int index = 0;
        			int everage = balls / colNum;
        
        			while (balls > 0){
        				if(colTotalSum[index] <= everage){
        					tmpRow[index] = colTotalSum[index];
        					balls -= colTotalSum[index];
                            if(balls==0)break;
        					index++;
							everage = balls / (colNum - index);
        				}else{
							tmpRow[index] = everage;
        					balls -= everage;
                            if(balls==0)break;
        					index++;
        					if(colNum != index)everage = balls / (colNum - index);
        				}
        			}
        
        			double chisq_tmp;
        			size_t df=0;
        			vector<double> p_null(colNum, 1.0/colNum);
        			ChisquareTest1DNoPValue(tmpRow, p_null, colNum, chisq_tmp, df);
        			boundHeteroChisq += chisq_tmp;
        			tmpRow.clear();
        		}
        		//cout<<boundHeteroChisq - m_colSumChisq[k] <<"  "<<m_totalObservedChisq<<endl;
        		result = (boundHeteroChisq - m_colSumChisq[k] +1e-07 > m_totalObservedChisq) ?
        		"to skip entire branch" : "not-to-skip";
			}
        
        
        
        
					 if(result == "to skip entire branch") {
						 // cout << "Suppose to skip";
						 m_skip = true;
						 // result = "not-to-skip";
        
						 /*
						  cout << "skip: actual chisq="
						  << m_boundChisqs[0] + m_boundChisqs[1] - m_nullHomoChisq
						  << " bound=" << boundHeteroChisq
						  << " observed=" << m_observedHeteroChisq
						  << endl;
						  */
					 }
        }
        return result;
    }
    //////
    //Hua added, Apr 21 2014
	void EMTFunctionalChisq::setRowAndColSum(vector<vector<int> > &row, vector<vector<int> > &col, const vector<TransitionTable> & Cs){
        m_requiredRowSums.resize(row.size());
        m_requiredColSums.resize(col.size());
//m_nullFisherProb.resize(Cs.size());
        
        for (size_t k=0; k<row.size(); k++) {
            m_requiredRowSums[k].resize(row[k].size());
            m_requiredColSums[k].resize(col[k].size());
            for (size_t i=0; i<row[k].size(); i++) {
                m_requiredRowSums[k][i] = row[k][i];
            }
            for (size_t i=0; i<col[k].size(); i++) {
                m_requiredColSums[k][i] = col[k][i];
            }
        }
    }


	double exact_functional_test(const vector< vector<int> > & C,
                             const string & discrepancy_measure)
{
	TransitionTable tt;
    tt.setTransitionTable(C);

	vector<TransitionTable> Cs(1, tt);

    EMTEnumerator e;
    
    EMTFunctionalChisq v;
    double P = exact_multi_table_test(Cs, v, e);
    if(v.getExtremeness()==true){
        P = 1 - P;
    }
    return P;
}
