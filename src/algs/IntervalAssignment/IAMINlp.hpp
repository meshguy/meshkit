// IAMINnl.hpp
// Interval Assignment for Meshkit
//
// ipopt mixed-integer solution
// The idea is the optimal solution will be an integer one, if one exists
// Define some region around the relaxed solution and search for an integer solution
//

#ifndef MESHKIT_IA_IAMINLP_HP
#define MESHKIT_IA_IAMINLP_HP

class IAData;
class IPData;
class IASolution;
class IANlp;

#include "MKVersion.h"
#include "IpTNLP.hpp"
using namespace Ipopt;

class IAMINlp : public TNLP
{
public:
  /** default constructor */
  IAMINlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr); 

  /** default destructor */
  virtual ~IAMINlp();

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
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);
  //@}


private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IAMINlp();
  IAMINlp(const IAMINlp&);
  IAMINlp& operator=(const IAMINlp&);
  //@}
  
  // input data
  const IAData *data;
  const IPData *ip_data;
  // solution data
  IASolution *solution;
  
  const bool debugging;
  const bool verbose; // verbose debugging
  
  // internally used methods
  // contributions of one variable to the objective function and gradient
  // underlying function 


  // r functions: if x>I then x-I / I else I-x / x
  Number eval_r_i(const Number& I_i, const Number& x_i); 
  Number eval_grad_r_i(const Number& I_i, const Number& x_i); 
  Number eval_hess_r_i(const Number& I_i, const Number& x_i); 

  // s functions: r weighted by x: r*x
  Number eval_s_i(const Number& I_i, const Number& x_i); 
  Number eval_grad_s_i(const Number& I_i, const Number& x_i); 
  Number eval_hess_s_i(const Number& I_i, const Number& x_i); 

  // Capital functions, l-p norms of the lowercase functions
public:
  Number eval_R_i(const Number& I_i, const Number& x_i); 
private:
  Number eval_grad_R_i(const Number& I_i, const Number& x_i); 
  Number eval_hess_R_i(const Number& I_i, const Number& x_i); 

  Number eval_S_i(const Number& I_i, const Number& x_i); 
  Number eval_grad_S_i(const Number& I_i, const Number& x_i); 
  Number eval_hess_S_i(const Number& I_i, const Number& x_i); 
  

};

#endif
