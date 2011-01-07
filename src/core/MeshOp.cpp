#include <assert.h>

#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit 
{
    
MeshOp::MeshOp(const MeshOp &mesh_op) throw(Error) : mkCore(mesh_op.mk_core())
{
  for (MOVector::const_iterator vit = mesh_op.parents().begin(); vit != mesh_op.parents().end(); vit++)
    (*vit)->add_child(this);
  
  for (MOVector::const_iterator vit = mesh_op.children().begin(); vit != mesh_op.children().end(); vit++)
    add_child(*vit);
}
  
MeshOp::MeshOp(MKCore *mkcore) throw(Error) : mkCore(mkcore) 
{
  if (mkcore) {
    mkcore->add_child(this);
    add_child(mkcore);
  }
}

MeshOp::~MeshOp() throw(Error)
{
  if (!opParents.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<MeshOp*>::reverse_iterator vit = opParents.rbegin(); vit != opParents.rend(); vit++)
      (*vit)->remove_child(this);
  }
  
  if (!opChildren.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<MeshOp*>::reverse_iterator vit = opChildren.rbegin(); vit != opChildren.rend(); vit++)
      remove_child(*vit);
  }

  assert(opParents.size() == 1 && opParents[0] == mkCore &&
         opChildren.size() == 1 && opChildren[0] == mkCore);
}

void MeshOp::add_parent(MeshOp *par) throw(Error)
{
  assert(("Either initializing mkCore" && this == mkCore && par == this) ||
         ("Or par isn't already a parent" &&
          std::find(opParents.begin(), opParents.end(), par) == opParents.end()));
  
    // if currently no parents, I'll no longer be a root, so remove from mkCore's children
  if (opParents.size() == 1 && opParents[0] == mkCore) {
    mkCore->remove_child(this);
    opParents[0] = par;
  }
  else
    opParents.push_back(par); 
}

void MeshOp::remove_parent(MeshOp *par) throw(Error)
{
  std::vector<MeshOp*>::iterator vit = std::find(opParents.begin(), opParents.end(), par);
  if (vit != opParents.end()) opParents.erase(vit);
  else throw Error(MK_FAILURE, "Parent operation not in parent list.");

    // if no more parents, make a child of the mkCore
  if (opParents.empty()) mkCore->add_child(this);
}

void MeshOp::add_child(MeshOp *child) throw(Error)
{
  assert(("Either initializing mkCore" && this == mkCore && child == this) ||
         ("Or child isn't already a child" &&
          std::find(opChildren.begin(), opChildren.end(), child) == opChildren.end()));

    // if no children, will no longer be a leaf, so remove as parent of mkCore
  if (opChildren.size() == 1 && opChildren[0] == mkCore) {
    mkCore->remove_parent(this);
    opChildren[0] = child;
  }
  else
    opChildren.push_back(child);

    // add myself as a parent of child
  child->add_parent(this);
}

void MeshOp::remove_child(MeshOp *child) throw(Error)
{
  std::vector<MeshOp*>::iterator vit = std::find(opChildren.begin(), opChildren.end(), child);
  if (vit != opChildren.end()) opChildren.erase(vit);
  else throw Error(MK_FAILURE, "Child operation not in child list.");

    // if no more children, I'm now a leaf, add as parent of mkCore
  if (opChildren.empty()) {
    mkCore->add_parent(this);
    opChildren.push_back(mkCore);
  }

    // remove myself as parent of child
  child->remove_parent(this);
}

  //! Setup function, called in reverse order before execute
void MeshOp::setup() throw(Error) 
{
  for (std::vector<MeshOp*>::iterator vit = opParents.begin(); vit != opParents.end(); vit++)
    (*vit)->setup();
    
  setup_this();
}
  
  //! Execute function, called in forward order after setup
void MeshOp::execute() throw(Error) 
{
  execute_this();

  for (std::vector<MeshOp*>::iterator vit = opChildren.begin(); vit != opChildren.end(); vit++)
    (*vit)->execute();
}

}

  
