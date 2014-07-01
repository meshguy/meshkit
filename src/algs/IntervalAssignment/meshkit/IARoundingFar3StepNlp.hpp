// IARoundingFarNnl.hpp
// Interval Assignment for Meshkit
//
// ipopt mixed-integer solution
// The idea is the optimal solution will be an integer one, if one exists
// Define some region around the relaxed solution and search for an integer solution
//

#ifndef MESHKIT_IA_IAROUNDINGFARP3STEPNLP_HP
#define MESHKIT_IA_IAROUNDINGFARP3STEPNLP_HP

class IAData;
class IPData;
class IASolution;
#include "meshkit/IANlp.hpp"
#include "meshkit/IAWeights.hpp"

#include "IpTNLP.hpp"

class IARoundingFar3StepNlp : public TNLP
{
public:
  /** default constructor */
  IARoundingFar3StepNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr); 

  /** default destructor */
  virtual ~IARoundingFar3StepNlp();
  
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
  IARoundingFar3StepNlp();
  IARoundingFar3StepNlp(const IARoundingFar3StepNlp&);
  IARoundingFar3StepNlp& operator=(const IARoundingFar3StepNlp&);
  //@}
  
  // input data
  const IAData *data;
  const IPData *ipData;
  IAWeights h0; // h0, as weights in IARoundingNLP
  IAWeights hp; // h+
  IAWeights hm; // h-
  // model variables indices
  const int x01_start;
  const int xp_start;
  const int xm_start;
  const int sum_even_start;
  const int x_constraint_start;
  const int problem_n, problem_m, base_n, base_m; // Index?

  // f' is discontinuous at integer values
  // Two choices: 
  // ZERO. give f'' as constant, zero, which is correct but hides the discontinuities
  // ROUNDED. give f'' as an average value, which might help the solver ..
  enum HessOptions {ZERO, ROUNDED};
  const HessOptions hess_option;

  // solution data
  IASolution *solution;
  
  
  // implemented using an overlay over an IANlp
  IANlp baseNlp;  
  
  const bool debugging;
  const bool verbose; // verbose debugging

  // utility
  double f_x_value( double I_i, double x_i ) const;
  double get_f_xl(int i) const; // obj function value at xl
  double get_f_xh(int i) const; // obj function value at xh
};

#endif
