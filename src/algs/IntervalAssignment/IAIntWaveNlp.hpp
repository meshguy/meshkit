// IAIntWaveNlp.hpp
// Interval Assignment for Meshkit
//
// ipopt mixed-integer solution
// The idea is the optimal solution will be an integer one, if one exists
// Define some region around the relaxed solution and search for an integer solution
//
// use some non-linear function as a constraint, so that the only feasible solutions are integer

#ifndef MESHKIT_IA_IAINTHESSNLP_HP
#define MESHKIT_IA_IAINTHESSNLP_HP

#include "IANlp.hpp"

// from Ipopt
#include "IpTNLP.hpp"

#include <vector>
#include <map>

namespace MeshKit 
{
    
class IAData;
class IPData;
class IASolution;

class IAIntWaveNlp : public TNLP
{
  // first set of functions required by TNLP
public:
  /** default constructor */
  IAIntWaveNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr,
               const bool set_silent = true); 

  /** default destructor */
  virtual ~IAIntWaveNlp();

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

// extra stuff not required by TNLP
  
private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IAIntWaveNlp();
  IAIntWaveNlp(const IAIntWaveNlp&);
  IAIntWaveNlp& operator=(const IAIntWaveNlp&);
  //@}
  
protected:
  // input data
  const IAData *data;
  const IPData *ipData;

  // solution data
  IASolution *solution;
  
  
  // implemented using an overlay over an IANlp
  IANlp baseNlp;  
  
  const int base_n, base_m;
  const int problem_n, problem_m;
  const int wave_even_constraint_start;
  const int wave_int_constraint_start;
  
  double f_x_value( double I_i, double x_i );

  const bool silent;
  const bool debugging;
  const bool verbose; // verbose debugging
  
  struct SparseMatrixEntry
	{
	  // position in matrix
	  // column must be less than row, j <= i
	  int i; // row
	  int j; // col
	  
	  // order in sequence for hessian values array, etc.
	  int k;
	
	  static int n; // matrix is n x n
	  int key() const { return i * n + j; }
	  
	  SparseMatrixEntry(const int iset, const int jset, const int kset);
	  SparseMatrixEntry() : i(-1), j(-1), k(-1) {} // bad values if unspecified
	};
	
	typedef std::map<int, SparseMatrixEntry> SparseMatrixMap; 
	
	SparseMatrixMap hessian_map;  // sorted by key
	std::vector< SparseMatrixEntry > hessian_vector; // from 0..k

	void add_hessian_entry( int i, int j, int &k );
	void build_hessian();
	int get_hessian_k( int i, int j );
  void print_hessian(); // debug
  
  // derived classes must define these
  // the first two should have a max of 1 for integer/even values, and be less than that elsewhere
  virtual double eval_g_int_x( const double x ) = 0;
  virtual double eval_g_int_s( const double s ) = 0;
  virtual double eval_jac_int_x( const double x ) = 0;
  virtual double eval_jac_int_s( const double s ) = 0;
  virtual double eval_hess_int_x( const double x ) = 0;
  virtual double eval_hess_int_s( const double s ) = 0;

};

} // namespace MeshKit

#endif
