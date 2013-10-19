#include "LPSolveClass.hpp"
#include "lp_lib.h"
#include <algorithm>

namespace MeshKit
{
//---------------------------------------------------------------------------//
// construction function for LPSolveClass class
LPSolveClass::LPSolveClass()
{
	

}



//---------------------------------------------------------------------------//
// deconstruction function for LPSolveClass class
LPSolveClass::~LPSolveClass()
{

}

//setup the objective function  : cx
void LPSolveClass::SetupObj(vector<double> left, double const_value)
{
	obj_const = const_value;
	coefficients.resize(left.size());
	for (unsigned int i = 0; i < left.size(); i++)
		coefficients[i] = left[i];
}

//setup the equal constraints   : Ax = b
void LPSolveClass::SetupEqu(vector<vector<double> > left, vector<double> right){

	assert(left.size()==right.size());
	A_equ.resize(left.size());
	b_equ.resize(right.size());
	for (unsigned int i = 0; i < left.size(); i++){
		A_equ[i].resize(left[i].size());
		b_equ[i] = right[i];
		for (unsigned int j = 0; j < left[i].size(); j++)
			A_equ[i][j] = left[i][j];
	}
}

//setup the inequal constraints: always, Ax <= b
void LPSolveClass::SetupInEqu(vector<vector<double> > left, vector<double> right){

	assert(left.size()==right.size());
	A_inequ.resize(left.size());
	b_inequ.resize(right.size());
	for (unsigned int i = 0; i < left.size(); i++){
		A_inequ[i].resize(left[i].size());
		b_inequ[i] = right[i];
		for (unsigned int j = 0; j < left[i].size(); j++)
			A_inequ[i][j] = left[i][j];
	}
}
//setup the const constraint
void LPSolveClass::SetupConst(vector<int> right)
{
	b_const.resize(right.size());
	for (unsigned int i = 0; i < right.size(); i++)
		b_const[i] =right[i];
}

//Execute function
int LPSolveClass::Execute()
{
	/*
	std::cout << "---------------------------------\n";
	std::cout << "objective function\n";
	for (unsigned int i = 0; i < coefficients.size(); i++)
		std::cout << coefficients[i] << "\t";
	std::cout << "\nConstant Value = " << obj_const << std::endl;

	std::cout << "---------------------------------\n";
	std::cout << "Equality Constraints\n";	
	for (unsigned int i = 0; i < A_equ.size(); i++){
		//std::cout << "Row index = " << i << "\t\t";
		for (unsigned int j = 0; j < A_equ[i].size(); j++)
			std::cout << A_equ[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "b\n";
	for (unsigned int i = 0; i < b_equ.size(); i++)
		std::cout << b_equ[i] << "\t";
	std::cout << "\n";


	std::cout << "---------------------------------\n";
	std::cout << "InEquality Constraints\n";	
	for (unsigned int i = 0; i < A_inequ.size(); i++){
		//std::cout << "Row index = " << i << "\t\t";
		for (unsigned int j = 0; j < A_inequ[i].size(); j++)
			std::cout << A_inequ[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "b\n";
	for (unsigned int i = 0; i < b_inequ.size(); i++)
		std::cout << b_inequ[i] << "\t";
	std::cout << "\n";
	*/

	lprec *lp;
	int Ncol = coefficients.size(), *colno = NULL, j, ret = 0;
	REAL *row = NULL;
	
	/* We will build the model row by row
     So we start with creating a model with 0 rows and n columns */

	lp = make_lp(0, Ncol);
	if (lp == NULL)
		ret = 1;/* couldn't construct a new model... */
		
	if (ret == 0){
		//let us name our variables
		std::string s = "x";
		for (int i = 0; i < Ncol; i++){
			std::stringstream out;
			out << i;
			s = s + out.str();
			char *cpy = new char[s.size()+1] ;
			strcpy(cpy, s.c_str());			
			set_col_name(lp, i+1, cpy);
		}

		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
    	row = (REAL *) malloc(Ncol * sizeof(*row));
		if ((colno == NULL) || (row == NULL))
      		ret = 2;
	}

	set_add_rowmode(lp, TRUE);
	//add the equation constraints
	if (ret == 0){
		/* makes building the model faster if it is done rows by row */
		if (A_equ.size() > 0){
			for (unsigned int i = 0; i < A_equ.size(); i++){//loop over the rows of equality constraints
				for (unsigned int j = 0; j < A_equ[i].size(); j++){//loop over the columns of equality constraints
					colno[j] = j+1;//add the j-th column to lpsolve
					row[j] = A_equ[i][j];
				}
				/* add the row to lpsolve */
				if(!add_constraintex(lp, A_equ[i].size(), row, colno, EQ, b_equ[i]))
					ret = 2;
			}
		}
	}
	
	//add the inequality constraints
	if (ret == 0){
		/* makes building the model faster if it is done rows by row */
		if (A_inequ.size() > 0){
			for (unsigned int i = 0; i < A_inequ.size(); i++){//loop over the rows of inequality constraints
				for (unsigned int j = 0; j < A_inequ[i].size(); j++){//loop over the columns of inequality constraints
					colno[j] = j+1;//add the j-th column to lpsolve
					row[j] = A_inequ[i][j];
				}
				/* add the row to lpsolve */
				if(!add_constraintex(lp, A_inequ[i].size(), row, colno, LE, b_inequ[i]))
					ret = 3;
			}
		}
	}

	//add the const constraint	
	if (ret == 0){
		if (b_const.size()>0){
			for (unsigned int i = 0; i < b_const.size(); i++){
				if (b_const[i] > 0){
					for (unsigned int j = 0; j < b_const.size(); j++){
						if (i == j){
							colno[j] = j+1;//add the j-th column to lpsolve
							row[j] = 1.0;						
						}				
						else{
							colno[j] = j+1;//add the j-th column to lpsolve
							row[j] = 0.0;
						}
					}
					if(!add_constraintex(lp, b_const.size(), row, colno, EQ, b_const[i]))
						ret = 4;		
				}
			}
		}
	}

	//set the variables to be integer
	if (ret == 0){
		for (int i = 0; i < Ncol; i++)
			set_int(lp, i+1, TRUE);
	}
	
	/* rowmode should be turned off again when done building the model */
	set_add_rowmode(lp, FALSE);	
	//add the objective function
	if (ret == 0){
		//set the objective function
		for (unsigned int i = 0; i < coefficients.size(); i++){
			colno[i] = i+1;
			row[i] = coefficients[i];
		}
		//set the objective in lpsolve
		if(!set_obj_fnex(lp, coefficients.size(), row, colno))
      		ret = 4;
	}

	//set the objective to minimize
	if (ret == 0){
		set_minim(lp);

		/* just out of curioucity, now show the model in lp format on screen */
    	/* this only works if this is a console application. If not, use write_lp and a filename */
    	write_LP(lp, stdout);

		/* I only want to see important messages on screen while solving */
    	set_verbose(lp, IMPORTANT);

    	/* Now let lpsolve calculate a solution */
    	ret = solve(lp);
    	if(ret == OPTIMAL)
      		ret = 0;
    	else
      		ret = 5;
	}

	//get some results
	if (ret == 0){
		/* a solution is calculated, now lets get some results */

    	/* objective value */
    	std::cout << "Objective value: " << get_objective(lp) << std::endl;

		/* variable values */
    	get_variables(lp, row);

		/* variable values */
		variables.resize(Ncol);
		for(j = 0; j < Ncol; j++)
			variables[j] = row[j];

		/* we are done now */
	}
	else{
		std::cout << "The optimal value can't be solved for linear programming, please check the constraints!!\n";
		exit(1);

	}
		
	
	std::cout << "print the result\t # of line segments is \n";
	for (int i = 0; i < Ncol; i++)
		std::cout << "index = " << i << "\t# = " << variables[i] << std::endl;

	/* free allocated memory */
  	if(row != NULL)
    	free(row);
  	if(colno != NULL)
    	free(colno);

	/* clean up such that all used memory by lpsolve is freed */
	if (lp != NULL)
		delete_lp(lp);

	return ret;
}

void LPSolveClass::GetVariables(vector<int> &var){
	var.resize(variables.size());
	for(unsigned int i = 0; i < variables.size(); i++)
		var[i] = variables[i];

}

}


