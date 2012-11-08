// IAInterface.hpp stubbed version

#ifndef MESHKIT_IA_INTERFACE_HP
#define MESHKIT_IA_INTERFACE_HP

/* stubbed includes
#include "meshkit/IAVariable.hpp"
#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/Interface.hpp"
*/

#include "IAVariable.hpp" // the one from stubs, also defines ModelEnt

#include <set>
#include <vector>

namespace MeshKit {

class IAData;
class IASolution;

/** \class IAInterface IAInterface.hpp "meshkit/IAInterface.hpp"
 * \brief The class used in MeshKit for Interval Assignment.
 *
 * Instances of this class are tools. The problem is to set up and solving the number of mesh edges 
 * to place on model entities and entity features. When solved, each curve can be meshed 
 * independently, and treated as fixed when meshing each surface or volume containing it. 
 * Different mesh schemes have different requirements (constraints), such as a mapped surface
 * needs opposite sides to have equal numbers of mesh edges (intervals).
 * The number of mesh edges on one curve is a variable; there may be additional variables.
 * The goal (objective function) is to have mesh edges close to the user-desired sizes.
 *
 * Construction of IAInterface does not cause the construction of variables.
 * Destruction of IAInterface does destroy the underlying variables.
 * IAInterface owns an IAVariable, for creation and deletion.
 * But the ModelEnt can request a variable for itself, and lets the interface know when 
 * it is no longer wanted.
 * ModelEnt keeps handles (pointers) to the variables it cares about.
 * \nosubgrouping
 */
 
class IAInterface //stubbed : public MeshScheme // register it with SchemeFactory
{
public:
   /** \name Constructor/destructor
     */
    /**@{*/

    /** \brief Constructor; model entity can be missing, in which case it's retrieved or created
     *
     * \param MKCore instance
     * \param MEntVector 
     */
// stubbed  IAInterface(MKCore *mkcore, const MEntVector &me_vec = MEntVector()) : MeshScheme(mkcore, me_vec){}

      /** \brief Destructor, destroys IAVariables
     */
  virtual ~IAInterface();
      /**@}*/
      
  /** \name Set Up
     */
    /**@{*/

   /** \name Variables
     */
    /**@{*/
   
     /** \brief Create a variable.
     * Variables are created before constraints.
     * It is OK if a variable is not in any constraint, but constraints are defined by their variables.
     * \param ModelEnt* model entity: this variable corresponds to the number of intervals
     *  on that model entity. If NULL, the variable corresponds to anything else that has
     *  meaning to the caller, and the caller will have to keep track of the variable 
     *  in order to access its solution value later.
     */
  IAVariable *get_variable( ModelEnt* model_entity = NULL, bool create_if_missing = true );
     /** \brief Create a variable and assign it a firmness and goal
     * \param IAVariable::Firmness The required fidelity of the solution to the goal. 
     * If HARD, then it is required that the solution equals the goal; goal should be integer.
     * If SOFT, usual case, try to get close to the goal.
     * If LIMP, we don't care how far the solution is from the goal.
     * \param double goal The desired number of intervals for this variable.
     * The goal may be non-integer, but we assume the solution must be a natural number, 
     * i.e. an integer >= 1.
     */
  // there is a reason we don't provide default parameters here, don't combine with above version.
  IAVariable *create_variable( ModelEnt* model_entity, IAVariable::Firmness set_firmness, double goal_value);
  
    /** \brief Get const_iterators over the variables. 
	*/
  typedef std::vector< IAVariable* > VariableVec;
  VariableVec::const_iterator variables_begin() const {return variables.begin();}
  VariableVec::const_iterator variables_end() const {return variables.end();}
  
     /** \brief Destroy a variable. If a variable is not explicitly destroyed, it will be
     * destroyed on IAInterface tool destruction.
     */
  void destroy_variable( IAVariable* ia_variable );
      /**@}*/

     /** \name Constraints
     */
    /**@{*/

     /** \brief Containers for variables for specifying constraints.
     */  
  //typedef std::vector<ModelEnt*> MEVec;
  // Types.hpp defines std::vector<ModelEnt*> MeshKit::MEntVector
  typedef std::vector<IAVariable*> IAVariableVec;
     /** \brief Convert container of ModelEnts to a container of IAVariables.
     * MEVec is an indirect way of specifying the model entities's variables.
     */   
    // convert vector of ModelEntities into vector of IAVariables
  IAVariableVec make_constraint_group( const MEntVector &model_entity_vec );
  
      /** \brief Constrain that the sum of the number of intervals 
      * on one side is equal to the number on the other side. E.g. when mapping a surface,
      * the opposite sides require equal intervals.
      */
  void constrain_sum_equal( const IAVariableVec &side_one, const IAVariableVec &side_two );
      /** \brief Constrain that the sum of the number of intervals is an even number, 
      * i.e. 2k for some integer k. E.g. for the curves bounding an unstructured quad mesh. 
      */
  void constrain_sum_even( const IAVariableVec &sum_even_vars );
      /** \brief More constraint types may be implemented here.
      */
  //... additional constraint types...
     /**@}*/
    /**@}*/

      /** \brief Main function that graph calls. Inherited from MeshScheme. 
      */
  virtual void setup_this();

  /** \name Solve the problem
     */
    /**@{*/

      /** \brief Main function that graph calls. Inherited from MeshScheme. 
      * find solution satisfying all the constraints
      * assign the solution to the variables
      * May be unsuccessful if the problem is over-specified. 
      * (how can callers detect failure? exception throw? no valid solution value in variables?) 
      */
  virtual void execute_this();
    /**@}*/

  /**\brief Get class name */
  static const char* name() 
    { return "IntervalAssignment"; }

  /**\brief Function returning whether this scheme can mesh entities of 
   *        the specified dimension.
   *\param dim entity dimension
   */
// stubbed
//  static bool can_mesh(iBase_EntityType dim)
//    { return iBase_VERTEX <= dim && iBase_REGION >= dim; }

  /** \brief Function returning whether this scheme can mesh the specified entity
   * 
   * Used by MeshOpFactory to find scheme for an entity.
   * \param model_ent ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
// stubbed
//  static bool can_mesh(ModelEnt *model_ent)
//      { return can_mesh((iBase_EntityType)model_ent->dimension()); }
    
  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
// stubbed
//  static const moab::EntityType* output_types();

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
// stubbed
//  virtual const moab::EntityType* mesh_types_arr() const
//    { return output_types(); }

  /** \brief Print the problem that was defined
   */
  void print_problem() const;
  
  /** \brief Destroy all the variables and constraints
   */
  void destroy_data();
  
private:
  /** \brief Internal representation of the data specifying the Interval Assignment problem.
   */
  // data
  VariableVec variables;
  typedef std::vector< VariableVec > VariableVecVec;
  VariableVecVec sumEqualConstraints1, sumEqualConstraints2; // one for each side
  VariableVecVec sumEvenConstraints;
  // ... additional types of constraints ...
  
  /** \brief Find the global interface index of the given variable
  */
  int variable_to_index(const IAVariable* var) const;
  /** \brief Find the ind'th variable.
  */
  IAVariable *index_to_variable(int ind) const;


  /** \name Subdivide into independent sub-problems
     */
    /**@{*/

    /** \brief Subdivide the problem into independent subproblems, in the sense that the
    * solution to one subproblem is not affected in any way by the solution to another.
    * I.e. independent if have no constraints or variables in common. 
    * A given subproblem is complete, in that it contains all the variables for each of 
    * its constraints, and also all the constraints for each of its variables. 
    */

  /** \name Represent a sub-set for a sub-problem 
     */  
  typedef std::set<int> IndexSet;
  /** \name Represent a sub-set for a sub-problem, vector for variable or constraint
     */  
  typedef std::vector< IndexSet > IndexSetVec;
  typedef std::vector< int > IndexVec;
  typedef std::vector< IndexVec > IndexVecVec;

  /** \brief Set = [0,1,...k-1]
  */
  void make_set_0_to_nm1( IndexSet &index_set, const size_t k); 
  /** \brief Vector = [0,1,...k-1]
  */
  void make_vec_0_to_nm1( IndexVec &index_vec, const size_t k); 
  /** \brief Vector = [-1,-1,...-1] of length k
  */
  void make_vec_unset( IndexVec &index_vec, const size_t k );

  /** \brief Represent the dependency of the problem using indicator sets.
  */
  struct VariableConstraintDependencies
  {
    /** brief 
    * IndexVecVec constraint_variables the ith constraint_variables are the indices of the variables in the ith constraints.
    * IndexVecVec variable_constraints the ith variable_constraints are the indices of the constraints containing the ith variable.
    */
    IndexVecVec constraintVariables, variableConstraints;
    void print() const; // debug
  };  

  /** \brief Build a representation of the dependency of the problem using indicator sets.
  * \param constraintVariables (member) Output. 
  * \param variableConstraints (member) Output. 
  */
  void set_variable_constraint_indices(VariableConstraintDependencies &var_con_dep);

  /** \brief Underlying workhorse to build a representation of the dependency of the 
  *   problem using indicator sets, for one constraint.
  * \param constraintVariables (member) Output. 
  * \param variableConstraints (member) Output. 
  * \param int i_start: Input. The following vector of constraints is indexed starting at i_start.
  * \param VariableVecVec: Input. Constraints, each entry is a vector specifying one constraint.
  */
  void set_variable_constraint_indices( VariableConstraintDependencies &var_con_dep,
  										const int i_start, 
                                        const VariableVecVec &variable_vec_vec );

  /** \brief Collections of indicator sets of variables and constraints that 
  * define a subproblem, or define the global problem, or define the remaining part of 
  * the global problem
  */
  class ProblemSets
  {
   public:
    IndexSet constraintSet, variableSet, hardVariableSet;
    IndexVec varMap; // map from local indices to indices in other set, i.e sub->global or global->sub
    bool empty() const {return constraintSet.empty() && variableSet.empty();}
    // default constructor makes empty sets, with an empty map, suitable for sub_problem initialization
    
    /** \brief Create a subset by chasing variables or constraints. Top of recursion - call this version. 
    * "this" is the global set.
    */
    void find_dependent_set( const VariableConstraintDependencies &var_con_dep, ProblemSets &subsets );
 
                                       
    /** \brief 
    * This is the recursive implementation; call "find_dependent_set" instead.
    * Add the variable (index) to the sub-problem, mark it as removed from the larger
    * problem, and recursively build the sub-problem by chasing its dependent constraints.
    * "this" is the global set.
    */
    void find_variable_dependent_set( const int variable_j, 
                                      const VariableConstraintDependencies &var_con_dep, 
                                      ProblemSets &subsets );
  
    /** \brief 
    * This is the recursive implementation; call "find_dependent_set" instead.
    * Add the constraint (index) to the sub-problem, mark it as removed from the larger
    * problem, and recursively build the sub-problem by chasing its variables.
    * "this" is the global set.
    */
    void find_constraint_dependent_set( const int constraint_i, 
                                        const VariableConstraintDependencies &var_con_dep, 
                                        ProblemSets &subsets );

    void print() const; // debug
  };
  /** \brief Initialize the global set defining the global problem.
  */
  void make_global_set(ProblemSets &problem_sets);

  /** \brief Convert a globally-indexed constraint into a locally-indexed constraint
  */
  void global_to_sub_side( const VariableVec &global_constraint, IndexVec &global_var_map, 
                           IndexVec &local_constraint, int &rhs ) const;

  /** \brief Represent a subproblem, its IAData and its mapping back to the global problem
  */
  struct SubProblem
  {
    static int max_id;
    int id;
    IAData *data;
    ProblemSets problemSets;
    IASolution *solution;
    void print() const; // debug
    SubProblem();
    virtual ~SubProblem();
  };
  typedef std::vector<SubProblem*> SubProblemVec;
  
  /** \brief Convert set-based definition of problem into an IAData based representation.
  */
  void fill_problem( ProblemSets &sub_sets, SubProblem *sub_problem, IndexVec &global_var_map ) const;

  /** \brief Build independent subproblems
  */
  void subdivide_problem(SubProblemVec &subproblems);
  void subdivide_problem_one(std::vector<IAData*> &subproblems); // just make one problem
  
  /**@}*/

  /** \brief For a subproblem, find a solution for the number of intervals for each variable.
  */
  bool solve_subproblem( SubProblem *subproblem );
  
  /** \brief Assign the found solution to the IAVariables.
  */
  void assign_solution( SubProblem *subproblem );
  
  static const bool debugging;
  
};

} // namespace MeshKit

#endif

/*
todo:
It shouldn't be a singleton, since there are cases I can envision where you'd want a local one.

Since IAInterface is a GraphNode, it will get destructed when the graph goes away.

- tim

On 10/26/2012 05:38 PM, Mitchell, Scott A wrote:

I'm a little confused about mkCore constructing an IAInterface.
How, if at all, is it destructed?

Should there be just one IAInterface, or many?

- Scott
*/

