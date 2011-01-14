#include <assert.h>

#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit 
{
    
MeshOp::MeshOp(const MeshOp &mesh_op) 
  : GraphNode(mesh_op), mkCore(mesh_op.mk_core())
{
}
  
MeshOp::MeshOp(MKCore *mkcore, const MEVector &me_vec) : mkCore(mkcore) 
{
  if (mkcore) {
    add_child(mkcore);
  }

  if (!me_vec.empty()) {
    for (MEVector::const_iterator vit = me_vec.begin(); vit != me_vec.end(); vit++)
      moab::Range &rit = meSelection[*vit];
  }
}

bool MeshOp::add_modelent(ModelEnt *model_ent) 
{
  MESelection::iterator sit = meSelection.find(model_ent);
  if (sit != meSelection.end()) return false;
  
  moab::Range &rtmp = meSelection[model_ent];
  return true;
}

bool MeshOp::remove_modelent(ModelEnt *model_ent) 
{
  MESelection::iterator sit = meSelection.find(model_ent);
  if (sit != meSelection.end()) {
    meSelection.erase(sit);
    return true;
  }
  else return false;
}

MeshOp::~MeshOp()
{
}

}

  
