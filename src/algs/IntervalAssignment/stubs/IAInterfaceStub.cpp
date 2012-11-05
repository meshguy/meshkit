// IAInterfaceStub.cpp stubbed version, without most of meshkit
// Interval Assignment Data for Meshkit
// Interface to the rest of Meshkit

#include "IAInterface.hpp" // the one in stubs = this dir
#include "IASolver.hpp"
#include "IADataBuilder.hpp"
#include "IASolution.hpp"

// stubbed #include "meshkit/ModelEnt.hpp" // from MeshKit

#include <vector>
#include <iterator>
#include <set>
#include <cstdio>
#include <assert.h>
#include <math.h>

namespace MeshKit 
{
  
int ModelEnt::max_stub_id = 0;   // stubbed
  
const bool IAInterface::debugging = false;
  
//static IAVariable counter
unsigned int IAVariable::numVariables(0);

// stubbed  
//moab::EntityType IAInterface_types[] = {moab::MBMAXTYPE};
//const moab::EntityType* IAInterface::output_types() 
//  { return IAInterface_types; }
    
void IAInterface::setup_this()
{
  ; //nothing for now
}
    
IAVariable *IAInterface::get_variable( ModelEnt* model_entity, bool create_if_missing )
{
  // already exists?
  if (model_entity && model_entity->ia_var())
  {
    // return without setting values
    return model_entity->ia_var();
  }
  
  if (model_entity && create_if_missing) 
  {
    // default values
    double goal(1.);
    Firmness firm(DEFAULT);
    
    const Firmness me_firm( model_entity->interval_firmness() );
    const double me_edge_size( model_entity->mesh_interval_size() );
    if ( me_firm == DEFAULT || me_edge_size <= 0.) 
    {
      // ModelEnt was unset, so use default values
      ;
    }
    else
    {
      const double me_size = model_entity->measure();
      assert(me_edge_size > 0.);
      goal = me_size / me_edge_size;
      firm = model_entity->interval_firmness(); 
    }
    return create_variable( model_entity, firm, goal);  
  }
  return NULL;
}


IAVariable *IAInterface::create_variable( ModelEnt* model_entity, IAVariable::Firmness set_firm, double goal_value)
{
  IAVariable *ia_var = NULL;

  // already exists?
  if (model_entity && model_entity->ia_var())
  {
    // set new values and return
    ia_var = model_entity->ia_var();
    ia_var->set_firmness(set_firm);
    ia_var->set_goal(goal_value);
    return ia_var;
  }
      
  ia_var = new IAVariable(model_entity, set_firm, goal_value);
  assert( ia_var->uniqueId + 1 == IAVariable::numVariables );
  variables.push_back(ia_var);
  assert( IAVariable::numVariables == variables.size() );

  return ia_var;
}

IAInterface::IAVariableVec IAInterface::make_constraint_group( const MEntVector &model_entity_vec )
{
  IAVariableVec result( model_entity_vec.size() );
  for (unsigned int i = 0; i < model_entity_vec.size(); ++i)
  {
    assert(model_entity_vec[i]);
    IAVariable *v = get_variable( model_entity_vec[i], true );
    result[i] = v;
  }
  return result;
}

void IAInterface::destroy_variable( IAVariable* ia_variable )
{
  if (!ia_variable)
    return;
  // model_entity shouldn't point to ia_variable anymore
  ModelEnt *me = ia_variable->get_model_ent();
  if (me && me->ia_var() == ia_variable)
    me->ia_var(NULL);
  
  variables[ ia_variable->uniqueId ] = NULL;
  delete ia_variable;    
}


void IAInterface::destroy_data()
{
  // destroy remaining variables
  for (unsigned int i = 0; i < variables.size(); ++i)
  {
    if ( variables[i] )
    {
      assert( variables[i]->uniqueId == i);
      destroy_variable( variables[i] );
    }
  }
  variables.clear();

  // reset global variable id? Yes if this is a singleton class, which it is right now.
  IAVariable::numVariables = 0;

  // clear constraints
  sumEqualConstraints1.clear();
  sumEqualConstraints2.clear();
  sumEvenConstraints.clear();
  // ... additional types of constraints ...
}

IAInterface::~IAInterface()
{
  destroy_data();
}

void IAInterface::constrain_sum_even( const IAVariableVec &sum_even_vars )
{
  sumEvenConstraints.push_back( sum_even_vars ); // vector copy
}

void IAInterface::constrain_sum_equal( const IAVariableVec &side_one, const IAVariableVec &side_two )
{
  sumEqualConstraints1.push_back(side_one); //vector copy of side_one
  sumEqualConstraints2.push_back(side_two); //vector copy of side_two
}


// for global vector of variables
int IAInterface::variable_to_index(const IAVariable* var) const
{
  return var->uniqueId;
}

IAVariable *IAInterface::index_to_variable(int ind) const
{
  assert(ind < (int) variables.size());
  return variables[ind];
}

void IAInterface::make_set_0_to_nm1( IndexSet &index_set, const size_t k )
{
  std::pair< IndexSet::iterator, bool> prior =  index_set.insert( 0 );
  assert(prior.second);
  for ( size_t i = 0; i < k; ++i)
  {
    prior.first = index_set.insert( prior.first, (int) i );
    // assert(prior.second); // second = insert-success is not assigned by this version of insert!
  }
}

void IAInterface::make_vec_0_to_nm1( IndexVec &index_vec, const size_t k )
{
  index_vec.clear();
  index_vec.reserve(k);
  for (size_t i = 0; i < k; ++i)
    index_vec.push_back( (int) i);
}

void IAInterface::make_vec_unset( IndexVec &index_vec, const size_t k )
{
  index_vec.clear();
  index_vec.reserve(k);
  for (int i = 0; i < (int) k; ++i)
    index_vec.push_back(-1);
}

void IAInterface::set_variable_constraint_indices(
            VariableConstraintDependencies &var_con_dep, 
            const int i_start, 
            const VariableVecVec &variable_vec_vec )
{
  unsigned int j;
  int i;
  for (j = 0, i = i_start; j < variable_vec_vec.size(); ++j, ++i)
  {
    for (unsigned int k = 0; k < variable_vec_vec[j].size(); ++k )
    {
      const IAVariable *var = variable_vec_vec[j][k];
      const int v = variable_to_index( var ); 
      var_con_dep.constraintVariables[i].push_back( v );        
      assert( v < (int) var_con_dep.variableConstraints.size() );
      var_con_dep.variableConstraints[v].push_back( i );
    }
  }
}
   
void IAInterface::set_variable_constraint_indices(VariableConstraintDependencies &var_con_dep)
{
  // create vectors of the length of constraints and variables, but containing of empty vectors
  var_con_dep.constraintVariables.resize( sumEqualConstraints1.size() + sumEvenConstraints.size() );
  var_con_dep.variableConstraints.resize( variables.size() );

  assert( sumEqualConstraints1.size() == sumEqualConstraints2.size() );
  const int sum_even_start = (int) sumEqualConstraints1.size();
  //zzyk hardsets
  // build vector map of variable-to-constraints
  // build vector map of constraint-to-variables
  set_variable_constraint_indices( var_con_dep, 0, sumEqualConstraints1 );
  set_variable_constraint_indices( var_con_dep, 0, sumEqualConstraints2 );
  set_variable_constraint_indices( var_con_dep, sum_even_start, sumEvenConstraints );
}

void IAInterface::ProblemSets::find_variable_dependent_set( 
         const int variable_j, 
         const VariableConstraintDependencies &var_con_dep,
         ProblemSets &subsets )
{
  // return if we've already added this variable to the subset = removed it from this set
  IndexSet::iterator j = variableSet.find(variable_j);
  if ( j == variableSet.end() )
    return;
  
  // add the variable to the sub-problem
  subsets.variableSet.insert(*j);
  
  // remove the variable from the big problem
  variableSet.erase(j);
  
  // recursively find the dependent constraints
  IndexVec::const_iterator i;
  assert( variable_j < (int) var_con_dep.variableConstraints.size() ); 
  const IndexVec &v = var_con_dep.variableConstraints[variable_j];
  for ( i = v.begin(); i != v.end(); ++i )
    find_constraint_dependent_set( *i, var_con_dep, subsets );
}


void IAInterface::ProblemSets::find_constraint_dependent_set( 
         const int constraint_i, 
         const VariableConstraintDependencies &var_con_dep,
         ProblemSets &subsets )
{
  // return if we've already added this constraint to the subset = removed it from this set
  // todo, could speed this up using some sort of mark on the entity
  IndexSet::iterator i = constraintSet.find(constraint_i);
  if ( i == constraintSet.end() )
    return;
  
  // add the constraint to the sub-problem
  subsets.constraintSet.insert(*i);

  // remove the constraint from the global problem
  constraintSet.erase(i);

  // recursively find the dependent variables
  IndexVec::const_iterator j;
  assert( constraint_i < (int) var_con_dep.constraintVariables.size() );
  const IndexVec &v = var_con_dep.constraintVariables[constraint_i];
  for ( j = v.begin(); j != v.end(); ++j )
    find_variable_dependent_set( *j, var_con_dep, subsets );
}

void IAInterface::ProblemSets::find_dependent_set( 
  const VariableConstraintDependencies &var_con_dep,
  ProblemSets &subsets )
{
  if ( !constraintSet.empty() )
  {
    find_constraint_dependent_set( *constraintSet.begin(), var_con_dep, subsets );
  }
  else if ( !variableSet.empty() )
  {
    find_variable_dependent_set( *variableSet.begin(), var_con_dep, subsets );
  }
}
   
void IAInterface::make_global_set(ProblemSets &problem_sets)
{
  make_set_0_to_nm1( problem_sets.constraintSet, sumEqualConstraints1.size() + sumEvenConstraints.size() );
  make_set_0_to_nm1( problem_sets.variableSet, variables.size() ); 
  make_vec_unset( problem_sets.varMap, variables.size() );
  // todo hardsets, 
  // move variables from variableSet to hard_variableSet if the variable is hard-set
}
  
void IAInterface::global_to_sub_side( 
  const VariableVec &global_constraint, 
  IndexVec &global_var_map,
  IndexVec &local_constraint, int &rhs ) const
{
  for (unsigned int j = 0; j < global_constraint.size(); ++j )
  {
    IAVariable *v = global_constraint[j];
    if (v->get_firmness() == MeshKit::HARD )
      rhs -= floor( v->get_goal() + 0.5 ); // goal should be integer already for hardsets
    else
    {
      int local_index = global_var_map[ v->uniqueId ];
      assert( local_index >= 0 );
      // assert( local_index < variableSet.size() );
      local_constraint.push_back( local_index ); //variable_to_index( v, sub_variables ) );
    }
  }
}
  
void IAInterface::ProblemSets::print() const
{
  printf("ProblemSets %p:\n", this);
  if (!empty())
    printf("non-");
  printf("empty. %lu variables, %lu hardsets, %lu constraints\n", variableSet.size(), hardVariableSet.size(), constraintSet.size());
  printf("Variables: ");
  for ( IndexSet::const_iterator i = variableSet.begin(); i != variableSet.end(); ++i )
  {
    printf("%d ", *i);
  }
  printf("\n");
  printf("HardSet Variables: ");
  for ( IndexSet::const_iterator i = hardVariableSet.begin(); i != hardVariableSet.end(); ++i )
  {
    printf("%d ", *i);
  }
  printf("\n");  
  printf("Constraints: ");
  for ( IndexSet::const_iterator i = constraintSet.begin(); i != constraintSet.end(); ++i )
  {
    printf("%d ", *i);
  }
  printf("\n\n");
}

void IAInterface::VariableConstraintDependencies::print() const
{
  printf("VariableConstraintDependencies %p:\n",this);
  printf("%lu variables, %lu constraints\n", variableConstraints.size(), constraintVariables.size() );
  for (size_t i = 0; i < variableConstraints.size(); ++i)
  {
    printf("Variable %lu depends on constraints ", i);
    for (size_t j = 0; j < variableConstraints[i].size(); ++j)
    {
      printf("%d ", variableConstraints[i][j]);
    }
    printf("\n");
  }
  
  for (size_t i = 0; i < constraintVariables.size(); ++i)
  {
    printf("Constraint %lu depends on variables ", i);
    for (size_t j = 0; j < constraintVariables[i].size(); ++j)
    {
      printf("%d ", constraintVariables[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

  
void IAInterface::fill_problem( ProblemSets &sub_sets, 
  SubProblem *sub_problem, IndexVec &global_var_map ) const
{
  IADataBuilder data_builder( sub_problem->data ); 
  // number the variables from 0..n, by setting ia_index
  sub_problem->data->I.reserve(sub_sets.variableSet.size());
  sub_sets.varMap.reserve(sub_sets.variableSet.size());
  for (IndexSet::const_iterator i = sub_sets.variableSet.begin(); i != sub_sets.variableSet.end(); ++i)
  {
    size_t index = *i; // global index
    IAVariable *v = index_to_variable( (int) index );
    if ( v->get_firmness() != MeshKit::HARD )
    {
      // index maps
      global_var_map[index] = (int) sub_sets.varMap.size();
      sub_sets.varMap.push_back( (int) index );

      // add the goals
      data_builder.add_variable( v->goal );
    }
    else
    {
     // nothing to do, it is a constant, not a variable in the sub problem
    }
  }

  // convert constraints to IASolver format
  for (IndexSet::const_iterator i = sub_sets.constraintSet.begin(); i != sub_sets.constraintSet.end(); ++i )
  {
    unsigned int c = *i;
    if (c < sumEqualConstraints1.size()) // equal
    {
      IndexVec side_1, side_2;
      int rhs_1(0), rhs_2(0);
        
      global_to_sub_side( sumEqualConstraints1[c], global_var_map, side_1, rhs_1 );
      global_to_sub_side( sumEqualConstraints2[c], global_var_map, side_2, rhs_2 );
      
      data_builder.constrain_opposite_side_equal(side_1, side_2, rhs_1 - rhs_2 );   
      // todo check should be rhs_2 - rhs_1 ?
    }
    else // even 
    {
      c -= sumEqualConstraints1.size();
      assert( c >= 0 ); // not useful because it is unsigned
      assert( c < sumEvenConstraints.size() );

      IndexVec side;
      int rhs(0);
      global_to_sub_side( sumEvenConstraints[c], global_var_map, side, rhs ); // todo check if should be -rhs?
      data_builder.constrain_sum_even(side, rhs );   
    }
  }
}

void IAInterface::subdivide_problem(SubProblemVec &subproblems)
{
  // find independent subproblems, by index of constraint and variable
  ProblemSets global_prob_set;
  make_global_set(global_prob_set);
  if (debugging)
    global_prob_set.print(); 
  VariableConstraintDependencies var_con_dep;
  set_variable_constraint_indices(var_con_dep);
  if (debugging)
    var_con_dep.print(); 
  
  while ( !global_prob_set.empty() )
  {
    // define a sub-problem by indices
    SubProblem *subproblem = new SubProblem;
    subproblems.push_back( subproblem );
    global_prob_set.find_dependent_set(var_con_dep, subproblem->problemSets);
    
    if (debugging)
    {
      printf("subproblem peeled out:");
      subproblem->problemSets.print();
      printf("remaining global problem:");
      global_prob_set.print();
    }
    
    // Fill in solver with data from those indices    
    // convert the set into an IASolver, add data
    fill_problem( subproblem->problemSets, subproblem, global_prob_set.varMap );
    if (debugging)
      subproblem->print();
  }    
}

void IAInterface::subdivide_problem_one(std::vector<IAData*> &subproblems)
{
   // placeholder - just make one subproblem
  IAData *sub_problem = new IAData();
  subproblems.push_back(sub_problem);
  IADataBuilder data_builder( sub_problem ); // data  

  // this grabs the whole problem and converts it
  // variables numbered 0..n
  for(VariableVec::const_iterator i = variables.begin(); i != variables.end(); ++i)
  {
    IAVariable *v = *i;
    data_builder.add_variable( v->goal );
  }
  // equal constraints
  for (unsigned int i = 0; i < sumEqualConstraints1.size(); ++i)
  {
    std::vector<int> side_1, side_2;
    for (unsigned int j = 0; j < sumEqualConstraints1[i].size(); ++j)
    {
      IAVariable *v = sumEqualConstraints1[i][j];
      int index = variable_to_index(v);
      // sanity check that the indexing was correct
      assert( sub_problem->I[index] == v->goal ); 
      side_1.push_back(index);
    }
    for (unsigned int j = 0; j < sumEqualConstraints2[i].size(); ++j)
    {
      IAVariable *v = sumEqualConstraints2[i][j];
      int index = variable_to_index(v);
      side_2.push_back(index);
    }
    data_builder.constrain_opposite_side_equal(side_1, side_2, 0.);
  }
  // even constraints
  for (unsigned int i = 0; i < sumEvenConstraints.size(); ++i)
  {
    std::vector<int> side;
    for (unsigned int j = 0; j < sumEvenConstraints[i].size(); ++j)
    {
      IAVariable *v = sumEvenConstraints[i][j];
      int index = variable_to_index(v);
      side.push_back(index);
    }
    data_builder.constrain_sum_even(side, 0.);
  }

}

void IAInterface::print_problem() const
{
  // describe variables
  // describe constraints
}

void IAInterface::SubProblem::print() const
{
  IASolverTool solver_tool(data,solution);
  solver_tool.print();
  // describe variables
  // describe constraints
}

bool IAInterface::solve_subproblem( SubProblem *subproblem )
{
  IASolver solver( subproblem->data, subproblem->solution );
  return solver.solve();
}

void IAInterface::assign_solution( SubProblem *subproblem )
{
  // assign solution value from subproblem to model entities
  std::vector<double> &x_solution = subproblem->solution->x_solution; // shorthand
  for (unsigned int i = 0; i < x_solution.size(); ++i)
  {
    const double x = x_solution[i];
    // map index i to IAVariable, which will be different when we have non-full subproblems
    IAVariable *v = index_to_variable( subproblem->problemSets.varMap[i] );
    assert(v);
    assert( x > 0. );
    const int x_int = (int) floor( x + 0.5 );
    v->solution = x_int;
    if (v->get_model_ent())
    {
      v->get_model_ent()->mesh_intervals(x_int);
    }
  }
}

void IAInterface::execute_this()
{
  SubProblemVec subproblems;
  subdivide_problem(subproblems);
  for (unsigned int i = 0; i < subproblems.size(); ++i)
  {
    SubProblem *p = subproblems[i];
    if ( solve_subproblem(p) )
      assign_solution(p);
    else
    {
      ; // some error statement
    }
    delete p;
  }
  subproblems.clear();
  return;
}

  
IAInterface::SubProblem::SubProblem()
{
  data = new IAData;
  solution = new IASolution;
}
  
IAInterface::SubProblem::~SubProblem()
{
  delete data;
  data = NULL;
  delete solution;
  solution = NULL;
}

} // namespace MeshKit


