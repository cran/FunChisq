// AdjustPvalue.cpp 
//
// Curtis Luce, Joe Song
// Created: 2007
// Last modified: September 27, 2008
// Last modified: October 4, 2011.  File name changed from "adjusted_pvalue.cpp"
// Last modified: May 19, 2014. fix factorial function to calculate huge factorials.
// Last modified: Feb 4, 2015. Use boost multiprecision to calculate factorial.

//#include <iostream>
//#include <sstream>

#include <cmath>
#include "TransitionTable.h"
#include "boost/math/special_functions/factorials.hpp" //Hua added, Nov 10 2014
//#include "boost/multiprecision/cpp_dec_float.hpp" //Hua added, Feb 3 2015

using namespace boost::math;
//using namespace boost::multiprecision;

#if defined _WIN32 || defined _WIN64
#define isinf(x) (!_finite(x))
#endif

bool compare_abs(const int & a, const int & b){
  return abs(a)<abs(b);
}

double factorial_vector(vector<int> v){//Hua Add, May 2 2014
	//Hua added, Nov 10 2014, calculate factorial by boost
  sort(v.begin(), v.end(), compare_abs);
	double value=1.0;
	for (int i=0; i<v.size(); i++) {
		if(v[i]>0){
			value = value * factorial<double>(abs(v[i]));
		}else if(v[i]<0){
			value = value / factorial<double>(abs(v[i]));
		}else{
      
		}
  }
	return value;
	////
  
  //Hua added, Feb 3 2015, Use boost multiprecision int
//    cpp_dec_float_50 value = 1;
//    for (int i=0; i<(int)v.size(); i++) {
//        if(v[i]>0){
//            int x = abs(v[i]);
//            while (x>0) {
//                value = value * x;
//                x--;
//            }
//            //value = value * factorial<double>(abs(v[i]));
//        }else if(v[i]<0){
//            int x = abs(v[i]);
//            while (x>0) {
//                value = value / x;
//                x--;
//            }
//            //value = value / factorial<double>(abs(v[i]));
//        }else{//v[i]==0
//            
//        }
//    }
//    
//    return value.convert_to<double>();
    ////
}

//double factorial (double n)
//{
//	if(n==0)
//		return 1;
//    
//	double answer = 1;
//
//	for (int i = 1; i <= n; i++) {
//		answer *= i;
//        if(isinf(answer)) {
//            break;
//        }
//    }
//
//	return answer;
//}
