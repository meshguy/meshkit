#ifndef MESHOP
#define MESHOP

#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "iGeom.h"
#include "moab/Interface.hpp"
#include <vector>

namespace MeshKit {

class MKCore;
    

/** \class MeshOp MeshOp.hpp "meshkit/MeshOp.hpp"
 * \brief An operation that operates on mesh data
 *
 * The purpose of this class is to enable mesh-related operations to be added to
 * an operation graph.  Once the graph has been created, it can be executed by calling
 * the execute() member function.
 */
class MeshOp
{
public:

    //! Copy constructor
  MeshOp(const MeshOp &mesh_op) throw(Error);
  
    //! Bare constructor
  MeshOp(MKCore *mkcore) throw(Error);

    //! Destructor
  virtual ~MeshOp() throw(Error);
  
    //! Return parent meshops
  std::vector<MeshOp*> &parents();
      
    //! Return parent meshops
  const std::vector<MeshOp*> &parents() const;
      
    //! Return child meshops
  std::vector<MeshOp*> &children();
      
    //! Return child meshops
  const std::vector<MeshOp*> &children() const;
      
    //! Add parent meshop
  void add_parent(MeshOp *par) throw(Error);

    //! Remove parent meshop
  void remove_parent(MeshOp *par) throw(Error);

    //! Add child meshop
  void add_child(MeshOp *child) throw(Error);

    //! Remove child meshop
  void remove_child(MeshOp *child) throw(Error);

    //! Setup function, called in reverse order before execute
  virtual void setup() throw(Error);

    //! Execute function, called in forward order after setup
  virtual void execute() throw(Error);
  
    //! Pure-virtual setup function, must be implemented in children
  virtual void setup_this() throw(Error) = 0;

    //! Pure-virtual execute function, must be implemented in children
  virtual void execute_this() throw(Error) = 0;

    //! Get the associated MKCore object
  MKCore *mk_core() const;

protected:
    //! Parent MeshOp's
  std::vector<MeshOp*> opParents;

    //! Child MeshOp's
  std::vector<MeshOp*> opChildren;
  
    //! MKCore associated with this MeshOp
  MKCore *mkCore;
  
private:

};

inline MKCore *MeshOp::mk_core() const
{
  return mkCore;
}

inline std::vector<MeshOp*> &MeshOp::parents() 
{
  return opParents;
}

inline const std::vector<MeshOp*> &MeshOp::parents() const
{
  return opParents;
}

inline std::vector<MeshOp*> &MeshOp::children() 
{
  return opChildren;
}
    
inline const std::vector<MeshOp*> &MeshOp::children() const
{
  return opChildren;
}
    
}
#endif

  
