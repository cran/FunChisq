//
//  TransitionTableIO.cpp
//  gln
//
//  Created by Joe Song on 5/21/13.
//
//  Updated:
//    Feb 3, 2014. MS. Added two general member functions scan(ifstream)
//      and save(ofstream)

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TransitionTable.h"

bool TransitionTable::scan(const string & file)
{
  bool success = true;

	m_transitionTable.clear();

	ifstream ifs(file.c_str());
    
    if (ifs.is_open()) {
  		string line;
  		while(getline(ifs, line)){
  			istringstream iss(line);
  			int n;
  			std::vector<int> v;
  			while (iss >> n)
  			{
  				v.push_back(n);
  			}
  			m_transitionTable.push_back(v);
  		}
    } else {
        success = false;
    }
    
    return success;
}
