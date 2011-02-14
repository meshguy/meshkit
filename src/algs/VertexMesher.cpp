#include "meshkit/VertexMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "moab/EntityType.hpp"
#include "meshkit/iGeom.hh"
#include "meshkit/RegisterMeshOp.hpp"

namespace MeshKit 
{

// static registration of this  mesh scheme
static RegisterMeshOp<VertexMesher,true> INIT("VertexMesher", iBase_VERTEX, moab::MBVERTEX);
    
VertexMesher::VertexMesher(MKCore *mkcore, const MEntVector &me_vec) 
        : MeshScheme(mkcore, me_vec) 
{}

VertexMesher::~VertexMesher() 
{}

MeshOp *VertexMesher::factory(MKCore *mkcore, const MEntVector &me_vec) 
{
  if (mkcore->vertex_mesher()) return mkcore->vertex_mesher();
  else {
    VertexMesher *mesher = new VertexMesher(mkcore, me_vec);
    mkcore->vertex_mesher(mesher);
    return mesher;
  }
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

    //! Setup is a no-op, but must be provided since it's pure virtual
void VertexMesher::setup_this() 
{}
    
void VertexMesher::execute_this() 
{
  if (mentSelection.empty()) return;
  
    // generate vertices for each vertex
  unsigned int i;
  MEntSelection::iterator sit;
  for (sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
    double pos[3];
      // get the position
    (*sit).first->evaluate(pos[0], pos[0], pos[0], pos);
    moab::EntityHandle new_vert;
      // create the vertex
    moab::ErrorCode rval = mk_core()->moab_instance()->create_vertex(pos, new_vert);
    MBERRCHK(rval, "Trouble creating new vertex.");
      // add to meselection
    (*sit).second.insert(new_vert);
      // commit to ModelEnt
    (*sit).first->commit_mesh((*sit).second, COMPLETE_MESH);
  }
}

  
} // namespace MeshKit
