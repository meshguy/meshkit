#ifndef MESHKIT_NO_OP_HPP
#define MESHKIT_NO_OP_HPP

#include "meshkit/MeshOp.hpp"

namespace MeshKit {  

/**\brief Dummy Operation
 */
class NoOp : public MeshOp 
{
public:
  NoOp( MKCore* core ) : MeshOp( core ) {}
  void setup_this() { }
  void execute_this() { }
  virtual const moab::EntityType* mesh_types_arr() const
    { 
      static moab::EntityType e = moab::MBMAXTYPE;
      return &e;
    }
private:
  //!\brief no copying
  NoOp( const NoOp& );
  void operator=( const NoOp& );
};



} // namespace MeshKit

#endif // MESHKIT_NO_OP_HPP
