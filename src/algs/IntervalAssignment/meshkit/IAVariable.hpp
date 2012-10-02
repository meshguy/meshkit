// IAVariable.hpp
// Data object shared between Meshkit and Interval Assignment

#ifndef MESHKIT_IA_VARIABLE_HP
#define MESHKIT_IA_VARIABLE_HP

#include "meshkit/Types.hpp" // for Firmness

namespace MeshKit 
{

class ModelEnt;

/** \class IAVariable IAVariable.hpp "meshkit/IAVariable.hpp"
 * \brief Class defining a variable for interval assignment.
 * These are the data objects shared between Interval Assignment and the rest of MeshKit.
 *
 * These should be constructed/destructed by Meshkit through IAInterface, not directly.
 * \nosubgrouping
 */
 
class IAVariable
{
public:
  
  /** \brief Get/Set the model_entity directly associated with this variable, if any.
  * This makes the most sense for curves.
  * Schemes and Entities may have extra variables beyond or instead of this one, 
  * e.g. intervals for holes or sweep levels, "periodic" intervals for through holes,
  * boundary layers, skew control decomposition, bias.
  */
  void set_model_ent(ModelEnt* me) {modelEnt = me;}
  ModelEnt *get_model_ent() const {return modelEnt;}
  
  /** \brief Get/Set the firmness, the required fidelity of the solution to the goal.
  */
  // enum Firmness {DEFAULT, SOFT, HARD}; //defined already in Types.hpp
  typedef MeshKit::Firmness Firmness;
  
  void set_firmness(Firmness set_firm) { firmness = set_firm; }
  Firmness get_firmness() const {return firmness;}
  void set_hard() {firmness = HARD;}
  void set_soft() {firmness = SOFT;}
  void unset()    {firmness = DEFAULT;} 
  bool is_soft() const  {return firmness == SOFT;}
  bool is_hard() const  {return firmness == HARD;} 
  bool is_default() const {return firmness == DEFAULT;}
  
  /** \brief Get/Set the goal, the number of intervals we'd like this variable to have in
  * the solution.
  */
  void set_goal(double goal_value) {goal=goal_value;}
  double get_goal() const {return goal;}
  
  /** \brief Get the solution to the interval assignment problem for this variable.
  * The solution exists if IAInterface::execute_this has been called already and it was
  * successful, and the value returned will be >=1. 
  * Otherwise the value returned will be < 1.
  *
  * Only interval assignment should be setting the solution,
  * so we leave "set_solution" out of the public interface.
  */
  int get_solution() const {return solution;}
  
      /** \brief Constructor; mesh entity can be missing.
     * Call IAInterface::create_variable instead; see that for documentation.
     */
  //void IAVariable();
  IAVariable(ModelEnt *model_entity = NULL, Firmness set_firmness = SOFT, double goal_value = 0. );
  
private:
   /** \name Internal data representation.
     */
    /**@{*/
  Firmness firmness;
  /** \brief goal should be >=1., non-integer is OK.
  */
  double goal; 
  /** \brief solution <0 notsolved, >=1 solved
  */
  int solution; 
  /** \brief modelEnt, NULL is OK.
  */
  ModelEnt *modelEnt;
    /**@}*/
        
  /** \brief IAInterface is friend class to it can set the solution, and use its index.
  */
  static unsigned int numVariables;
  unsigned int uniqueId;
  friend class IAInterface;
};

inline
IAVariable::IAVariable(ModelEnt *model_entity, Firmness set_firmness, double goal_value )
  : firmness(set_firmness), goal(goal_value), solution(0), modelEnt(model_entity), uniqueId(numVariables++)
{}

} // namespace MeshKit
#endif
