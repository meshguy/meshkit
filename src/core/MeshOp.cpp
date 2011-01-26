#include <assert.h>

#include "meshkit/MeshOpFactory.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"

namespace MeshKit 
{
    
MeshOp::MeshOp(const MeshOp &mesh_op) 
        : GraphNode(mesh_op)
{
}
  
MeshOp::MeshOp(MKCore *mkcore, const MEVector &me_vec) 
        : GraphNode(mkcore)
{
  if (!me_vec.empty()) {
    for (MEVector::const_iterator vit = me_vec.begin(); vit != me_vec.end(); vit++)
      moab::Range &rit = meSelection[*vit];
  }
}

void MeshOp::setup_boundary() 
{
    // generic setup code; make sure there are MeshOp's for each of my bounding entities, and create
    // ones where necessary (based on default MeshOp for that dimension registered with MOF)

    // use separate MeshOps by dimension, since there may be multiple dimension entities in meSelection
  MeshOp *this_op[] = {NULL, NULL, NULL, NULL};

  for (MESelection::iterator mit = meSelection.begin(); mit != meSelection.end(); mit++) {
    ModelEnt *this_ent = (*mit).first;
    int dim = this_ent->dimension();
    if (0 == dim) continue;

    MEVector children;
    this_ent->get_adjacencies(dim-1, children);
    for (MEVector::iterator chit = children.begin(); chit != children.end(); chit++) {
      MOVector chops;
      (*chit)->get_meshops(chops);
      if (chops.empty()) {
          // no meshop, need to get one associated with this
        if (!this_op[dim-1]) {
          this_op[dim-1] = MeshOpFactory::instance()->construct_meshop(dim-1);
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

  
