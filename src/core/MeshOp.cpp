#include <assert.h>

#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"

namespace MeshKit 
{
    
MeshOp::MeshOp(const MeshOp &mesh_op) 
        : GraphNode(mesh_op)
{
}
  
MeshOp::MeshOp(MKCore *mkcore, const MEntVector &me_vec) 
        : GraphNode(mkcore)
{
  if (!me_vec.empty()) {
    for (MEntVector::const_iterator vit = me_vec.begin(); vit != me_vec.end(); vit++)
      moab::Range &rit = mentSelection[*vit];
  }
}

void MeshOp::setup_boundary() 
{
    // generic setup code; make sure there are MeshOp's for each of my bounding entities, and create
    // ones where necessary (based on default MeshOp for that dimension registered with MOF)

    // use separate MeshOps by dimension, since there may be multiple dimension entities in mentSelection
  MeshOp *this_op[] = {NULL, NULL, NULL, NULL};

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    ModelEnt *this_ent = (*mit).first;
    int dim = this_ent->dimension();
    if (0 == dim) continue;

    MEntVector children;
    this_ent->get_adjacencies(dim-1, children);
    for (MEntVector::iterator chit = children.begin(); chit != children.end(); chit++) {
      MOpVector chops;
      (*chit)->get_meshops(chops);
      if (chops.empty()) {
          // no meshop, need to get one associated with this
        if (!this_op[dim-1]) {
          this_op[dim-1] = mk_core()->construct_meshop(dim-1);
          if (!this_op[dim-1]) throw Error(MK_MESHOP_NOT_FOUND, "No default meshop for this dimension.");
        }
        this_op[dim-1]->add_modelent(*chit);
      }
    }
  }
  
  for (int dim = 0; dim <= 3; dim++) {
    if (this_op[dim])
      mk_core()->insert_node(this_op[dim], this);
  }
}

bool MeshOp::add_modelent(ModelEnt *model_ent) 
{
  MEntSelection::iterator sit = mentSelection.find(model_ent);
  if (sit != mentSelection.end()) return false;
  
  moab::Range &rtmp = mentSelection[model_ent];
  return true;
}

bool MeshOp::remove_modelent(ModelEnt *model_ent) 
{
  MEntSelection::iterator sit = mentSelection.find(model_ent);
  if (sit != mentSelection.end()) {
    mentSelection.erase(sit);
    return true;
  }
  else return false;
}

MeshOp::~MeshOp()
{
}

}

  
