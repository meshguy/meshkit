// IARoundingFarNnl.hpp
// Interval Assignment for Meshkit
//
// ipopt mixed-integer solution
// The idea is the optimal solution will be an integer one, if one exists
// Define some region around the relaxed solution and search for an integer solution
//

#ifndef MESHKIT_IA_IAROUNDINGFARNLP_HP
#define MESHKIT_IA_IAROUNDINGFARNLP_HP

class IAData;
class IPData;
class IASolution;
#include "IANlp.hpp"

#include "IpTNLP.hpp"
using namespace Ipopt;

class IARoundingFarNlp : public TNLP
{
public:
  /** default constructor */
  IARoundingFarNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr); 

  /** default destructor */
  virtual ~IARoundingFarNlp();

  // extra for IA, not for TNLP
  bool randomize_weights_of_non_int();
  // return true if any weights were changed
  
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
  IARoundingFarNlp();
  IARoundingFarNlp(const IARoundingFarNlp&);
  IARoundingFarNlp& operator=(const IARoundingFarNlp&);
  //@}
  
  // input data
  const IAData *data;
  const IPData *ipData;
  // model
  std::vector<double> h0; // h0, as weights in IARoundingNLP
  std::vector<double> hp; // h+
  std::vector<double> hm; // h-
  
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
