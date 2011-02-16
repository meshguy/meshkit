#include "meshkit/VertexMesher.hpp"
#include "meshkit/MKCore.hpp"
#include "moab/EntityType.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/RegisterMeshOp.hpp"

namespace MeshKit 
{
    
 
// static registration of this  mesh scheme
moab::EntityType types[] = { moab::MBVERTEX, moab::MBMAXTYPE};
const moab::EntityType* VertexMesher::output_types() 
  { return types; }
    
VertexMesher::VertexMesher(MKCore *mkcore, const MEntVector &me_vec) 
        : MeshScheme(mkcore, me_vec) 
{}

VertexMesher::~VertexMesher() 
{
  mk_core()->vertex_mesher(NULL);
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
    MBERRCHK(rval, mk_core()->moab_instance());
      // add to meselection
    (*sit).second.insert(new_vert);
      // commit to ModelEnt
    (*sit).first->commit_mesh((*sit).second, COMPLETE_MESH);
  }
}

  
} // namespace MeshKit
