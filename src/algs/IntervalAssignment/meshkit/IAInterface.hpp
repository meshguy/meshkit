// IAVariable.hpp
// Define a variable for interval assignment
// Data object shared between Meshkit and Interval Assignment

#ifndef MESHKIT_IA_INTERFACE_HP
#define MESHKIT_IA_INTERFACE_HP

// #include "ModelEntity.hpp"
class ModelEntity;
class IASolver;
#include "IAVariable.hpp"
// class IAVariable;
#include <set>
#include <vector>

// to add an entity or variable to interval assignment 
// create an IAVariable object, with or without 


// IAInterface - owns an IAVariable, for creation and deletion.
// But the ModelEntity can request a variable for itself, and lets the interface know when it is no longer wanted.
// ModelEntity keeps handles (pointers) to the variables it cares about.
// most efficient interface: add variable as a constraint is added

class IAInterface // todo: derive from MeshScheme, register it with SchemeFactory
// todo: provide methods: setup_this, execute_this
{
public:
  // create a variable, without an associated constraint (yet)
  IAVariable *create_variable( ModelEntity* model_entity = NULL );
  // there is a reason we don't provide default parameters here
  IAVariable *create_variable( ModelEntity* model_entity, IAVariable::Firmness set_firmness, double goal_value);
  
  void destroy_variable( IAVariable* ia_variable );
  
  
  typedef std::vector<ModelEntity*> MEVec;
  typedef std::vector<IAVariable*> IAVariableVec;
  
    // convert vector of ModelEntities into vector of IAVariables
  IAVariableVec make_constraint_group( const MEVec &model_entity_vec ); //todo implement this
  
  void add_sum_even_constraint( const IAVariableVec &sum_even_vars );
  void add_sum_equal_constrataint( const IAVariableVec &side_one, const IAVariableVec &side_two );
  //... additional constraint types...

  // Main function that graph calls
  bool setup_this();

  // Main function that graph calls graph
  // find solution satisfying all the constraints
  // assign the solution to the variables
  // return true if successful
  bool execute_this();
  
  // constructor/destructor
  virtual ~IAInterface();
  IAInterface();
  
private:
  // data
  typedef std::set<IAVariable*> VariableSet;
  VariableSet variables;
  typedef std::vector< IAVariableVec > VariableSetVec;
  VariableSetVec sumEqualConstraints1, sumEqualConstraints2; // one for each side
  VariableSetVec sumEvenConstraints;
  // ... additional types ...
  
  // solution sub-methods
  // return the position of var in variables, from 0 .. n-1
  int variable_to_index(const IAVariable* var);
  typedef std::set<int> IndexSet;
  typedef std::vector< IndexSet > IndexSetVec;
  void make_0_to_nm1( const int k, IndexSet &index_set ); // populate index_set with 0..n-1
  
  
  // for subdivide_problem
  int variable_to_index(const IAVariable* var) const;
  void get_constraint_variable_indices( IndexSetVec &constraint_variables, 
                                       IndexSetVec &variable_constraints,
                                       const int i_start, 
                                       const VariableSetVec &variable_set_vec );
  void find_variable_dependent_set( const int variable_j, 
                                   IndexSetVec &constraint_variables,
                                   IndexSetVec &variable_constraints,
                                   IndexSet &constraint_set, IndexSet &variable_set, 
                                   IndexSet &sub_constraint_set, IndexSet &sub_variable_set);
  void find_constraint_dependent_set( const int constraint_i, 
                                     IndexSetVec &constraint_variables,
                                     IndexSetVec &variable_constraints,
                                     IndexSet &constraint_set, IndexSet &variable_set, 
                                     IndexSet &sub_constraint_set, IndexSet &sub_variable_set);
  
  
  void subdivide_problem(std::vector<IASolver*> &subproblems);
  bool solve_subproblem( IASolver *subproblem );
  void assign_solution( IASolver *subproblem );
  void destroy_subproblems(std::vector<IASolver*> &subproblems);
  
};

#endif
