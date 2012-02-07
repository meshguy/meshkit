#include <assert.h>

#include "meshkit/MeshOp.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"

//#include "meshkit/FBiGeom.hpp"

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
      mentSelection[*vit];
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
      if ((*chit)->is_meshops_list_empty()) {
          // no meshop, need to get one associated with this
        if (!this_op[dim-1]) {
          this_op[dim-1] = mk_core()->construct_meshop(dim-1);
          if (!this_op[dim-1]) throw Error(MK_MESHOP_NOT_FOUND, "No default meshop for this dimension.");
        }
        this_op[dim-1]->add_modelent(*chit);
      }
    }
  }

  // VertexMesher should be inserted after root node to be done before other meshops
  // It is because the VertexMesher instance is only one
  if (this_op[0]) {
    mk_core()->insert_node(this_op[0], this, mk_core()->root_node());
  }

  for (int dim = 1; dim <= 3; dim++) {
    if (this_op[dim])
      mk_core()->insert_node(this_op[dim], this);
  }
}

bool MeshOp::canmesh_vertex(ModelEnt *model_ent) 
{
  return (model_ent->dimension() == 0);
}

bool MeshOp::canmesh_edge(ModelEnt *model_ent)
{
  return (model_ent->dimension() == 1);
}

bool MeshOp::canmesh_face(ModelEnt *model_ent)
{
  return (model_ent->dimension() == 2);
}

bool MeshOp::canmesh_region(ModelEnt *model_ent)
{
  return (model_ent->dimension() == 3);
}

bool MeshOp::add_modelent(ModelEnt *model_ent) 
{
  MEntSelection::iterator sit = mentSelection.find(model_ent);
  if (sit != mentSelection.end()) return false;
  
  mentSelection[model_ent];

  // add meshop back to model ent
  model_ent->add_meshop(this);
  
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

void MeshOp::mesh_types(std::vector<moab::EntityType> &mesh_types) 
{
  const moab::EntityType* types = mesh_types_arr();
  for (int i = 0; types[i] != moab::MBMAXTYPE; ++i)
    mesh_types.push_back(types[i]);
}

#ifdef HAVE_FBIGEOM
void MeshOp::create_model_ents_from_previous_ops()
{
  // model ents are empty so far...
  GraphNode * prevNode = other_node(in_arcs());// assumes one incoming edge!!

  MeshOp * prevOp = reinterpret_cast<MeshOp*> (prevNode);
  std::string nameOp = prevNode->get_name();

  if (NULL==prevOp)
    return;
  // if the previous op is MBSplitOp or MBGeomOp, we know what to do
  if ( !(nameOp == "MBGeomOp" || nameOp == "MBSplitOp"))
  {
    return; // do not process yet other OPs
  }


  MEntSelection & mentsel = prevOp->me_selection();
  if (mentsel.empty())
    return;
  ModelEnt * firstMe = (*(mentsel.begin()) ).first;
  moab::Range prevModelSet = mentsel[firstMe];

  if (prevModelSet.size()!=1)
   return;

  // this range has one set that can serve as basis for
  // mesh based geometry
  // the model ents have no model tag on the sets!!
  iGeom::EntitySetHandle rootSet = (iGeom::EntitySetHandle) prevModelSet[0];

  int latestGeomIndex = mk_core()->initialize_mesh_based_geometry(rootSet);

  MEntVector model_ents;
  // get the model entities of dimension 2, created with the latest FBiGeom!!!
  mk_core()->get_entities_by_dimension(2, model_ents, latestGeomIndex);

  // add it to the mentSelection!!
  if (model_ents.size()==0)
    return; // nothing to mesh .... could be an infinite loop maybe we should abort
  for (unsigned int i = 0; i<model_ents.size(); i++)
  {
    moab::Range rr = mentSelection[model_ents[i]];
    rr.clear(); // just to use it , to avoid some warnings
  }

  /*// redo the setup for the op, as we added more entities
  setup_this();
  // we will have to also execute operations before this ...
  mk_core()->setup_and_execute();// will relaunch the whole execute// some are marked as
  // executed, so no worries*/

}
#endif

}

  
