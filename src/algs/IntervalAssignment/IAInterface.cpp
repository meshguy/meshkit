// IAInterface.cpp
// Interval Assignment Data for Meshkit
// Interface to the rest of Meshkit

#include "meshkit/IAInterface.hpp"
#include "IASolver.hpp"
#include "meshkit/ModelEnt.hpp" // from MeshKit
#include <vector>
#include <cstdio>
#include <assert.h>

namespace MeshKit 
{

IAVariable *IAInterface::create_variable( ModelEnt* model_entity )
{
  // already exists?
  if (model_entity && model_entity->ia_var())
  {
    // return without setting values
    return model_entity->ia_var();
  }
 
  return create_variable( model_entity, IAVariable::SOFT, 0. );  
}

void IAInterface::destroy_variable( IAVariable* ia_variable )
{
  if (!ia_variable)
    return;
  // model_entity shouldn't point to ia_variable anymore
  if (ia_variable->get_model_entity() && ia_variable->get_model_entity()->ia_var() == ia_variable)
    ia_variable->get_model_entity()->ia_var(NULL);
  
  variables.erase( ia_variable );
  delete ia_variable;    
}

IAInterface::~IAInterface()
{
  // destroy remaining variables
  // in reverse order for efficiency
  while (!variables.empty())
  {
    destroy_variable( * variables.rbegin() );
  }

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
  // add to set, giving the end of the set as a hint of where it goes
  variables.insert( ia_var);

  return ia_var;
}


void IAInterface::make_0_to_nm1( IndexSet &index_set, const int k )
{
  for (int i = 0; i < k; ++i)
    index_set.insert( i);
}

int IAInterface::variable_to_index(const IAVariable* var) const
{
  VariableSet::const_iterator var_pos = variables.find( const_cast<IAVariable*>(var) );
  int v = std::distance(variables.begin(), var_pos);
  assert( v >= 0 );
  assert( v < (int)variables.size() ); // == size means it wasn't found
  return v;
}

void IAInterface::get_constraint_variable_indices( IndexSetVec &constraint_variables, 
                                                  IndexSetVec &variable_constraints,
                                                  const int i_start, 
                                                  const VariableSetVec &variable_set_vec )
{
  unsigned int j;
  int i;
  for (j = 0, i = i_start; j < variable_set_vec.size(); ++j, ++i)
  {
    for (VariableSet::const_iterator k = variable_set_vec[j].begin(); k != variable_set_vec[j].end(); ++k)
    {
      int v = variable_to_index(*k); 
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
    /*
  // find independent subproblems, by index of constraint and variable
  // put all the constraints together, for generality
  //typedef std::set<int> IndexSet;
  IndexSet constraint_set, variable_set;
  IndexSet sub_constraints, sub_variables;
  assert( sumEqualConstraints1.size() == sumEqualConstraints2.size() );
  const int sum_equal_start = sumEqualConstraints1.size();
  make_0_to_nm1( constraint_set, sumEqualConstraints1.size() + sumEvenConstraints.size() );
  make_0_to_nm1( variable_set, variables.size() );
  IndexSetVec constraint_variables(constraint_set.size()), variable_constraints( variable_set.size() );
  get_constraint_variable_indices( constraint_variables, variable_constraints, 0, sumEqualConstraints1 );
  get_constraint_variable_indices( constraint_variables, variable_constraints, 0, sumEqualConstraints2 );
  get_constraint_variable_indices( constraint_variables, variable_constraints, sum_equal_start, sumEvenConstraints );
  

  while ( !constraint_set.empty() || !variable_set.empty())
  {
    // define a sub-problem by indices
    IndexSet sub_constraint_set, sub_variable_set;
    if ( !constraint_set.empty() )
    {
      find_constraint_dependent_set( *constraint_set.begin(), constraint_variables, variable_constraints, 
                                    constraint_set, variable_set, sub_constraints, sub_variables);
    }
    else 
    {
      assert( !variable_set.empty() );
      find_variable_dependent_set( *variable_set.begin(), constraint_variables, variable_constraints, 
                                    constraint_set, variable_set, sub_constraints, sub_variables);
    }

    // create a new sub-prolem, with data from those indices
    IASolver *sub_problem = new IASolver();

    
  //    -- zzyk: makes sure the below makes sense given the above problem description
  //  need to go from indices back to ia variables...
    // remove ia_index
    
  // convert the set into an IASolver, add data
  
  // number the variables from 0..n, by setting ia_index
    int j;
    ConstraintSet::const_iterator i;
    for (j = 0, i = sub_variables.begin(); i != sub_variables.end(); ++i, ++j)
  {
    IAVariable *v = index_to_variable( *i );
    if ( i->get_firmness() != IAVariable::HARD )
    {
      // add the goals
      sub_problem->I.push_back( i->goal );
      assert( sub_problem->num_variables() == j+1 );
    }
  }
  
  
  // convert equality constraints to IASolver  
  for (ConstraintSet::const_iterator i = sub_sum_equal_constraints.begin(); i != sub_sum_equal_constraints.end() )
  {
    IAVariableVec *s1 = &sumEqualConstraints1[*i];
    IAVariableVec *s2 = &sumEqualConstraints2[*i];
    std::vector<int> side_1, side_2;
    int rhs(0);
    
    for (int j = 0; j < s1->size(); ++j )
    {
      IAVariable *v = s1[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs -= floor( v->get_goal() + 0.5 ); // goal should be integer already for hardsets
      else
        side_1.push_back( v->ia_index );
    }
    for (int j = 0; j < s2->size(); ++j )
    {
      IAVariable *v = s2[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs += v->get_goal();
      else
        side_2.push_back( v->ia_index );
    }
    sub_problem->constrain_opposite_side_equal(side_1, side_2, rhs );   
  }
  
  // convert even constraints to IASolver  
  for (ConstraintSet::const_iterator i = sub_even_constraints.begin(); i != sub_even_constraints.end() )
  {
    IAVariableVec *s = &sumEvenConstraints[i];
    std::vector<int> side;
    int rhs(0);
    
    for (int j = 0; j < s->size(); ++j )
    {
      IAVariable *v = s[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs -= floor( v->get_goal() + 0.5 ); // goal should be integer already for hardsets
      else
        side.push_back( v->ia_index );
    }
    sub_problem->constrain_sum_even(side, rhs );   
  }
  
  //todo: find any variables that didn't have any constraints, solve those by simply rounding to the nearest integer.
  // for example, curves in triangle scheme sides
  }

    */

  IASolver *sub_problem = new IASolver();
  VariableSet::iterator vit;
  for (vit = variables.begin(); vit != variables.end(); vit++) 
  {
    IAVariable *v = index_to_variable( *vit );
    if ( vit->get_firmness() != IAVariable::HARD )
    {
      // add the goals
      sub_problem->I.push_back( vit->goal );
    }
  }

  // convert equality constraints to IASolver  
  unsigned int i;
  for (i = 0; i < sumEqualConstraints1.size(); i++)
  {
    IAVariableVec *s1 = &sumEqualConstraints1[i];
    IAVariableVec *s2 = &sumEqualConstraints2[i];
    std::vector<int> side_1, side_2;
    int rhs(0);
    
    for (int j = 0; j < s1->size(); ++j )
    {
      IAVariable *v = s1[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs -= floor( v->get_goal() + 0.5 ); // goal should be integer already for hardsets
      else
        side_1.push_back( v->ia_index );
    }
    for (int j = 0; j < s2->size(); ++j )
    {
      IAVariable *v = s2[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs += v->get_goal();
      else
        side_2.push_back( v->ia_index );
    }
    sub_problem->constrain_opposite_side_equal(side_1, side_2, rhs );   
  }
  
  // convert even constraints to IASolver  
  for (ConstraintSet::const_iterator i = sub_even_constraints.begin(); i != sub_even_constraints.end() )
  {
    IAVariableVec *s = &sumEvenConstraints[i];
    std::vector<int> side;
    int rhs(0);
    
    for (int j = 0; j < s->size(); ++j )
    {
      IAVariable *v = s[j];
      if (v->get_firmness() == IAVariable::HARD )
        rhs -= floor( v->get_goal() + 0.5 ); // goal should be integer already for hardsets
      else
        side.push_back( v->ia_index );
    }
    sub_problem->constrain_sum_even(side, rhs );   
  }

  subproblems.push_back(sub_problem);
}

bool IAInterface::solve_subproblem( IASolver *subproblem )
{
/*
  for (int i=subproblems.size()-1; i>=0; --i)
  {
    delete subproblems[i];
  }
  subproblems.clear();
*/

  return false;
}

void IAInterface::assign_solution( IASolver *subproblem )
{
}

} // namespace MeshKit
