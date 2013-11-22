//-----------------------------------C++-------------------------------------//
// File: src/algs/Sweep/LPSolveClass.hpp
// Monday Feb 01 10:50 2012
// Brief: LPSolveClass class definition: do the LPSolveClass 
//        linear programming using lpsolve
//---------------------------------------------------------------------------//


#ifndef LPSOLVECLASS_HPP
#define LPSOLVECLASS_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <map>


using namespace std;


namespace MeshKit
{
//===========================================================================//
  /*!
   * \class LPSolveClass
   * \brief do the LPSolveClass to solve the linear programming problem
   * 
   * 
   */
//===========================================================================//

class LPSolveClass
{	
public:
    LPSolveClass();
    ~LPSolveClass();

	//setup the objective function and constraints
	void SetupObj(vector<double> left, double const_value = 0.0);
	void SetupInEqu(vector<vector<double> > left, vector<double> right);
	void SetupEqu(vector<vector<double> > left, vector<double> right);
	void SetupConst(vector<int> right);


    //Execute function
    int Execute();

	//return the variables
	void GetVariables(vector<int> &var);

   


private://private member functions
	

private://private member variable
	vector<double> coefficients;
	double obj_const;
	vector<vector<double> > A_inequ;
	vector<double> b_inequ;
	vector<vector<double> > A_equ;
	vector<double> b_equ;
	vector<int> b_const;

	vector<int> variables;
   
};

}
#endif
