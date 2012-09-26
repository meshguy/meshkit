// IAVariable.hpp
// Define a variable for interval assignment
// Data object shared between Meshkit and Interval Assignment

#ifndef MESHKIT_IA_INTERFACE_HP
#define MESHKIT_IA_INTERFACE_HP

// #include "ModelEnt.hpp"
#include "meshkit/IAVariable.hpp"
#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"

#include <set>
#include <vector>

// to add an entity or variable to interval assignment 
// create an IAVariable object, with or without 


// IAInterface - owns an IAVariable, for creation and deletion.
// But the ModelEnt can request a variable for itself, and lets the interface know when it is no longer wanted.
// ModelEnt keeps handles (pointers) to the variables it cares about.
// most efficient interface: add variable as a constraint is added

namespace MeshKit {

class ModelEnt;
class IASolver;

class IAInterface : public MeshScheme // register it with SchemeFactory
// todo: provide methods: setup_this, execute_this
{
public:
  // constructor/destructor
  virtual ~IAInterface();
  IAInterface(MKCore *mkcore, const MEntVector &me_vec = MEntVector());
  
  // create a variable, without an associated constraint (yet)
  IAVariable *create_variable( ModelEnt* model_entity = NULL );
  // there is a reason we don't provide default parameters here
  IAVariable *create_variable( ModelEnt* model_entity, IAVariable::Firmness set_firmness, double goal_value);
  
  void destroy_variable( IAVariable* ia_variable );
  
  
  typedef std::vector<ModelEnt*> MEVec;
  typedef std::vector<IAVariable*> IAVariableVec;
  
    // convert vector of ModelEntities into vector of IAVariables
  IAVariableVec make_constraint_group( const MEVec &model_entity_vec ); //todo implement this
  
//  void add_sum_even_constraint( const IAVariableVec &sum_even_vars );
//  void add_sum_equal_constraint( const IAVariableVec &side_one, const IAVariableVec &side_two );
  //... additional constraint types...

  // Main function that graph calls
  virtual void setup_this();

  // Main function that graph calls graph
  // find solution satisfying all the constraints
  // assign the solution to the variables
  // return true if successful
  virtual void execute_this();
  
  /**\brief Get class name */
  static const char* name() 
    { return "IntervalAssignment"; }

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
    { return iBase_VERTEX <= dim && iBase_REGION >= dim; }

  /** \brief Function returning whether this scheme can mesh the specified entity
   * 
   * Used by MeshOpFactory to find scheme for an entity.
   * \param model_ent ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *model_ent)
      { return can_mesh((iBase_EntityType)model_ent->dimension()); }
    
  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
  static const moab::EntityType* output_types();

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }

private:
  // data
  typedef std::set<IAVariable*> VariableSet;
  VariableSet variables;
  typedef std::vector< VariableSet > VariableSetVec;
  VariableSetVec sumEqualConstraints1, sumEqualConstraints2; // one for each side
  VariableSetVec sumEvenConstraints;
  // ... additional types ...
  
  // solution sub-methods
  // return the position of var in variables, from 0 .. n-1
  int variable_to_index(const IAVariable* var) const;
  IAVariable *index_to_variable(int ind) const;
  typedef std::set<int> IndexSet;
  typedef std::vector< IndexSet > IndexSetVec;
  void make_0_to_nm1( IndexSet &index_set, const int k); // populate index_set with 0..n-1
  
  
  // for subdivide_problem
  void get_constraint_variable_indices( IndexSetVec &constraint_variables, 
                                       IndexSetVec &variable_constraints,
                                       const int i_start, 
                                       const VariableSetVec &variable_set_vec );
  void find_variable_dependent_set( const int variable_j, 
                                    const IndexSetVec &constraint_variables,
                                    const IndexSetVec &variable_constraints,
                                    IndexSet &constraint_set, IndexSet &variable_set, 
                                    IndexSet &sub_constraint_set, IndexSet &sub_variable_set);
  
  void find_constraint_dependent_set( const int constraint_i, 
                                      const IndexSetVec &constraint_variables,
                                      const IndexSetVec &variable_constraints,
                                      IndexSet &constraint_set, IndexSet &variable_set, 
                                      IndexSet &sub_constraint_set, IndexSet &sub_variable_set);
  
  
  void subdivide_problem(std::vector<IASolver*> &subproblems);
  bool solve_subproblem( IASolver *subproblem );
  void assign_solution( IASolver *subproblem );
  void destroy_subproblems(std::vector<IASolver*> &subproblems);
  
};

} // namespace MeshKit

#endif
