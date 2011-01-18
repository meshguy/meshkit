#include "meshkit/VertexMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOpFactory.hpp"
#include "moab/EntityType.hpp"
#include "iGeom.hh"

namespace MeshKit 
{

// static registration of this  mesh scheme
int success = MeshOpFactory::instance()->register_meshop("VertexMesher", moab::MBVERTEX, VertexMesher::factory, NULL);
    
MeshOp *VertexMesher::factory(MKCore *mkcore, const MEVector &me_vec) 
{
  return new VertexMesher(mkcore, me_vec);
}

void VertexMesher::mesh_types(std::vector<moab::EntityType> &tps) 
{
  tps.push_back(moab::MBVERTEX);
}
    
bool VertexMesher::add_modelent(ModelEnt *model_ent) 
{
    // make sure this represents a geometric vertex
  if (0 != model_ent->dimension()) 
    throw Error(MK_WRONG_DIMENSION, "Wrong dimension entity added to VertexMesher.");

  return MeshOp::add_modelent(model_ent);
}

void VertexMesher::execute_this() 
{
  if (meSelection.empty()) return;
  
    // generate vertices for each vertex
  unsigned int i;
  MESelection::iterator sit;
  for (sit = meSelection.begin(); sit != meSelection.end(); sit++) {
    double pos[3];
      // get the position
    (*sit).first->closest(pos[0], pos[0], pos[0], pos);
    moab::EntityHandle new_vert;
      // create the vertex
    moab::ErrorCode rval = mkCore->moab_instance()->create_vertex(pos, new_vert);
    MBERRCHK(rval, "Trouble creating new vertex.");
      // add to meselection
    (*sit).second.insert(new_vert);
      // commit to ModelEnt
    (*sit).first->commit_mesh((*sit).second, COMPLETE_MESH);
  }
}

  
} // namespace MeshKit
