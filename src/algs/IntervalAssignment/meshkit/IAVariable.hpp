// IAVariable.hpp
// Define a variable for interval assignment
// Data object shared between Meshkit and Interval Assignment

#ifndef MESHKIT_IA_VARIABLE_HP
#define MESHKIT_IA_VARIABLE_HP

#include <cstring>

class ModelEntity;

class IAVariable
{
public:
  
  void set_model_entity(ModelEntity* me);
  ModelEntity *get_model_entity() const;
  
  enum Firmness {UNSET, SOFT, HARD};
  
  void set_firmness(Firmness set_firm) { firmness = set_firm; }
  Firmness get_firmness() const {return firmness;}
  void set_hard() {firmness = HARD;}
  void set_soft() {firmness = SOFT;}
  void unset()    {firmness = UNSET;} 
  bool is_soft() const  {return firmness == SOFT;}
  bool is_hard() const  {return firmness == HARD;} 
  bool is_unset() const {return firmness == UNSET;}
  
  void set_goal(double goal_value) {goal=goal_value;}
  double get_goal() const {return goal;}
  
  int get_solution() const {return solution;}
  // only interval assignment should be setting the solution
  // so we leave "set_solution" out of the public interface
  
  //void IAVariable();
  IAVariable(ModelEntity *model_entity = NULL, Firmness set_firmness = SOFT, double goal_value = 0. );
  
  
private:
  Firmness firmness;
  double goal; 
  int solution; // <0 notsolved, >=1 solved
  ModelEntity *modelEntity; // NULL is OK
  
  friend class IAInterface;
  int ia_index; // zzyk- remove this used during solution process, could replace with a temporary map is memory is an issue
};

inline
IAVariable::IAVariable(ModelEntity *model_entity, Firmness set_firmness, double goal_value )
: modelEntity(model_entity), goal(goal_value), firmness(set_firmness), solution(0)
{}

#endif
