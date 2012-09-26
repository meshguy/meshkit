// IANnl.hpp
// Interval Assignment for Meshkit
//
// This is the solver-translation from Meshkit to the ipopt library
// Here we provide the functions that ipopt needs to solve the optimization problem
//
// Adapted from: 
// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.hpp 1864 2010-12-22 19:21:02Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-09

#ifndef MESHKIT_IA_IANLP_HP
#define MESHKIT_IA_IANLP_HP

/** Optimization problem for getting a non-integer solution to the mapping and 
 *  other constraints that might cause a large change in interval sizes.
 * 
 * Solved using Ipopt
 *
 * Problem structure:
 *  minimize F = measure of (f_i)
            where sum x/I + I/x : x is assigned intervals, I is goal intervals
 *  st. g = sum side_1 - sum side_2 = 0 : for a mapping constraint
 *          sum side_1 + sum side_2 - sum side_3 >= 0 : triangle inequalities for midpoint subdivision, etc
 *          1 <= x <= infinity
 *
 * f has optimal value of x=I. The term I/x penalizes x<I, and x/I penalizes x>I.
 * f' = 1/I - I/x^2
 * f'' = 2I/x^3 , strictly positive
 *
 * F is some measure of the f_i, combo of 
 *  sum of the square of the f_i
 *  max(f_i)*number_of_f_i
 *  experiment with what works well
 *
 * g' = sum of non-zero coefficients of the constraints, since they are all linear.
 * g'' = 0
 *
 * Class hierarchy: The Meshkit layers are supposed to be independent, except they all 
 * depend on shared IAData. Implemented using a tower of classes, but multiple virtual 
 * inheritance would be closer to the conceptual model.
 *
 */

#include "MKVersion.h"

#include "IpTNLP.hpp"
using namespace Ipopt;

namespace MeshKit 
{
    
class IAData;
class IASolution;

class IANlp : public TNLP
{
public:
  /** default constructor */
  IANlp(const IAData *data_ptr, IASolution *solution_ptr); 

  /** default destructor */
  virtual ~IANlp();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x_init,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData *ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}
  
  int get_neleJac() const { return neleJac; }

private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IANlp();
  IANlp(const IANlp&);
  IANlp& operator=(const IANlp&);
  //@}
  
  // input data
  const IAData *data;
  static const int p_norm;
  // solution data
  IASolution *solution;
  int neleJac;
  
  const bool debugging;
  const bool verbose; // verbose debugging
  
  // internally used methods
  // contributions of one variable to the objective function and gradient
  // underlying function 


  // r functions: if x>I then x-I / I else I-x / x
  static Number eval_r_i(const Number& I_i, const Number& x_i); 
  static Number eval_grad_r_i(const Number& I_i, const Number& x_i); 
  static Number eval_hess_r_i(const Number& I_i, const Number& x_i); 

  // s functions: r weighted by x: r*x
  static Number eval_s_i(const Number& I_i, const Number& x_i); 
  static Number eval_grad_s_i(const Number& I_i, const Number& x_i); 
  static Number eval_hess_s_i(const Number& I_i, const Number& x_i); 

  // Capital functions, l-p norms of the lowercase functions
public:
  static Number eval_R_i(const Number& I_i, const Number& x_i); 
private:
  static Number eval_grad_R_i(const Number& I_i, const Number& x_i); 
  static Number eval_hess_R_i(const Number& I_i, const Number& x_i); 

  static Number eval_S_i(const Number& I_i, const Number& x_i); 
  static Number eval_grad_S_i(const Number& I_i, const Number& x_i); 
  static Number eval_hess_S_i(const Number& I_i, const Number& x_i); 
  

};

} // namespace MeshKit 

#endif
