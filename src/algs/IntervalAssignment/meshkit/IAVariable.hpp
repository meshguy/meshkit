// IAVariable.hpp
// Define a variable for interval assignment
// Data object shared between Meshkit and Interval Assignment

#ifndef MESHKIT_IA_VARIABLE_HP
#define MESHKIT_IA_VARIABLE_HP

#include <cstring>

namespace MeshKit 
{

class ModelEnt;

class IAVariable
{
public:
  
  void set_model_entity(ModelEnt* me);
  ModelEnt *get_model_entity() const;
  
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
  IAVariable(ModelEnt *model_entity = NULL, Firmness set_firmness = SOFT, double goal_value = 0. );
  
  
private:
  Firmness firmness;
  double goal; 
  int solution; // <0 notsolved, >=1 solved
  ModelEnt *modelEnt; // NULL is OK
  
  friend class IAInterface;
  int ia_index; // zzyk- remove this used during solution process, could replace with a temporary map is memory is an issue
};

inline
IAVariable::IAVariable(ModelEnt *model_entity, Firmness set_firmness, double goal_value )
  : firmness(set_firmness), goal(goal_value), solution(0), modelEnt(model_entity)
{}

} // namespace MeshKit
#endif
