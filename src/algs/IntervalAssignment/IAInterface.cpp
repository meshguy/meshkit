// IAInterface.cpp
// Interval Assignment Data for Meshkit
// Interface to the rest of Meshkit

#include "meshkit/IAInterface.hpp"
#include "IASolver.hpp"
#include "meshkit/ModelEnt.hpp" // from MeshKit
#include <vector>
#include <iterator>
#include <set>
#include <cstdio>
#include <assert.h>
#include <math.h>

namespace MeshKit 
{

//static IAVariable counter
unsigned int IAVariable::numVariables(0);
  
moab::EntityType IAInterface_types[] = {moab::MBMAXTYPE};
const moab::EntityType* IAInterface::output_types() 
  { return IAInterface_types; }
    
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

IAInterface::~IAInterface()
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

void IAInterface::make_0_to_nm1( IndexSet &index_set, const int k )
{
  for (int i = 0; i < k; ++i)
    index_set.insert( i);
}

int IAInterface::variable_to_index(const IAVariable* var) const
{
  return var->uniqueId;
}

IAVariable *IAInterface::index_to_variable(int ind) const
{
  assert(ind < (int) variables.size());
  return variables[ind];
}

void IAInterface::get_constraint_variable_indices( IndexSetVec &constraint_variables, 
                                                  IndexSetVec &variable_constraints,
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
      constraint_variables[i].insert( v );        
      assert( v < (int)variable_constraints.size() );
      variable_constraints[v].insert( i );
    }
  }
}

void IAInterface::find_variable_dependent_set( const int variable_j, 
                                              const IndexSetVec &constraint_variables,
                                              const IndexSetVec &variable_constraints,
                                              IndexSet &constraint_set, IndexSet &variable_set, 
                                              IndexSet &sub_constraint_set, IndexSet &sub_variable_set)
{
  // return if we've already added this variable
  IndexSet::iterator j = variable_set.find(variable_j);
  if ( j == variable_set.end() )
    return;
  
  // add the variable to the sub-problem
  sub_variable_set.insert(*j);
  
  // remove the variable from the big problem
  variable_set.erase(j);
  
  // recursively find the dependent constraints
  for (IndexSet::iterator i = variable_constraints[variable_j].begin(); i != variable_constraints[variable_j].end(); ++i)
    find_constraint_dependent_set( *i, constraint_variables, variable_constraints, 
                                  constraint_set, variable_set, sub_constraint_set, sub_variable_set);

}


void IAInterface::find_constraint_dependent_set( const int constraint_i, 
                                                const IndexSetVec &constraint_variables,
                                                const IndexSetVec &variable_constraints,
                                                IndexSet &constraint_set, IndexSet &variable_set, 
                                                IndexSet &sub_constraint_set, IndexSet &sub_variable_set)
{
  // return if we've already added this constraint
  IndexSet::iterator i = constraint_set.find(constraint_i);
  if ( i == constraint_set.end() )
    return;
  
  // add the constraint to the sub-problem
  sub_constraint_set.insert(*i);

  // remove the constraint from the big problem
  constraint_set.erase(i);

  // recursively find the dependent variables
  for (IndexSet::iterator j = constraint_variables[constraint_i].begin(); j != constraint_variables[constraint_i].end(); ++j)
    find_variable_dependent_set( *j, constraint_variables, variable_constraints, 
                                constraint_set, variable_set, sub_constraint_set, sub_variable_set);
}


void IAInterface::subdivide_problem(std::vector<IASolver*> &subproblems)
{

    /* todo: get subdivision working from src code in progress and copy it over */

   // placeholder - just make one subproblem
  IASolver *sub_problem = new IASolver();
  subproblems.push_back(sub_problem);

  // this grabs the whole problem and converts it
  // variables numbered 0..n
  for(VariableVec::const_iterator i = variables.begin(); i != variables.end(); ++i)
  {
    IAVariable *v = *i;
    sub_problem->I.push_back( v->goal );
  }
  // equal constraints
  for (unsigned int i = 0; i < sumEqualConstraints1.size(); ++i)
  {
    std::vector<int> side_1, side_2;
    for (unsigned int j = 0; j < sumEqualConstraints1[i].size(); ++j)
    {
      IAVariable *v = sumEqualConstraints1[i][j];
      int index = variable_to_index(v);
      side_1.push_back(index);
    }
    for (unsigned int j = 0; j < sumEqualConstraints2[i].size(); ++j)
    {
      IAVariable *v = sumEqualConstraints2[i][j];
      int index = variable_to_index(v);
      side_2.push_back(index);
    }
    sub_problem->constrain_opposite_side_equal(side_1, side_2, 0.);
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
    sub_problem->constrain_sum_even(side, 0.);
  }

}

bool IAInterface::solve_subproblem( IASolver *subproblem )
{
  return subproblem->solve();
}

void IAInterface::assign_solution( IASolver *subproblem )
{
  // assign solution value from subproblem to model entities
  for (unsigned int i = 0; i < subproblem->x_solution.size(); ++i)
  {
    const double x = subproblem->x_solution[i];
    // map index i to IAVariable, which will be different when we have non-full subproblems
    IAVariable *v = index_to_variable( i );
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
  std::vector<IASolver*> subproblems;
  subdivide_problem(subproblems);
  for (unsigned int i = 0; i < subproblems.size(); ++i)
  {
    IASolver* p = subproblems[i];
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

} // namespace MeshKit
