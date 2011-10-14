/*
 * MBSplitOp.cpp
 *
 *  Created on: Oct 2, 2011
 *      Author: iulian
 */

#include "meshkit/MBSplitOp.hpp"
#include "meshkit/ModelEnt.hpp"

#include "moab/Core.hpp"
#include "meshkit/FBiGeom.hpp"
#include "moab/FBEngine.hpp"

namespace MeshKit {

//Entity Type initialization for splitting; no mesh output
moab::EntityType MBSplitOp_tps[] = { moab::MBMAXTYPE }; // no mesh, really
const moab::EntityType* MBSplitOp::output_types()
{
  return MBSplitOp_tps;
}

MBSplitOp::MBSplitOp(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
{
}

MBSplitOp::~MBSplitOp()
{
}

//set up the crop/split of a mesh based surface
void MBSplitOp::set_options(int globalId, double dirx, double diry,
    double dirz, int closed)
{
  _globalId = globalId;
  /*for (int i=0; i<nPoints*3; i++)
   _polyline.push_back(polyline[i]);*/

  _direction[0] = dirx;
  _direction[1] = diry;
  _direction[2] = dirz;
  _closed = closed;
  return;
}
void MBSplitOp::add_points(double x, double y, double z)
{
  _polyline.push_back(x);
  _polyline.push_back(y);
  _polyline.push_back(z);
  return;
}
// model entities should be created on mesh based geometry
void MBSplitOp::setup_this()
{
  if (mentSelection.empty())
    return;
  // go through the map, to find the set with the given globalId
  // find the set of dimension 2, with the global id matching
  MEntSelection::iterator mit;

  moab::Interface * MBI = mk_core()->moab_instance();
  moab::EntityHandle mset;
  moab::Tag gid = mk_core()->moab_global_id_tag();
  moab::ErrorCode rval;
  int id = -1;
  for (mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    mset = (mit->first)->mesh_handle();
    // get the globalId tag

    rval = MBI->tag_get_data(gid, &mset, 1, &id);
    if (_globalId == id)
      break;
  }
  if (id != _globalId) {
    std::cout << " the face not found, abort\n";
    return;
  }
  //
  // get the one and only model entity, which should be of dimension 2, and geometrize it
  // this should be the model entity handle

  iGeom::EntityHandle gent = (mit->first)->geom_handle();
  if (gent != 0) {
    std::cout << "geometry exists on the model ent; Abort.\n";
    return; //
  }

  moab::FBEngine * pFacet = new moab::FBEngine(MBI,
      NULL, true);

  rval = pFacet->Init();
  MBERRCHK(rval, MBI);

  // the mset will be split according to the rule

  moab::EntityHandle newFace;
  rval = pFacet->split_surface_with_direction(mset, _polyline, _direction,
      newFace, _closed);

  MBERRCHK(rval, MBI);

  // at the end of this, moab database will contain new gsets, for the split eometry
  // the old pFacet should be cleaned and deleted
  //
  pFacet->delete_smooth_tags();
  delete pFacet;// will trigger a cleanup , including deleting the
  // smooth faces, edges, and tags;
  // now, the moab db would be ready for another round, including meshing...
  // maybe not all surfaces need to be meshed, only the "result" new face
  // it will get a new id too;

  // the result of the split operation will be a set encompassing all the
  // gsets that form the current topological model
  // and the newFace + g children of it
  // the root model could be used for a new GeomTopoTool, as in
  //  new  GeomTopoTool(mb, true, newRootModel);
  moab::EntityHandle newRootModel;
  rval = MBI->create_meshset(moab::MESHSET_SET, newRootModel);
  MBERRCHK(rval, MBI);
  // add to this root model set all the moabEntSets from model ents in this
  // MeshOp, and all children of the newFace !!!
  // collect in newRootModel all gsets...
  // they should form a "correct" topo model; (check!)
  for (mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
    mset = (mit->first)->mesh_handle();
    // get all children of this set
    moab::Range children;
    rval = MBI->get_child_meshsets(mset, children, 0); // all children
    MBERRCHK(rval, MBI);
    children.insert(mset); // insert the current face too
    rval = MBI->add_entities(newRootModel, children);
    MBERRCHK(rval, MBI);
  }
  // add the new face children too
  moab::Range children;
  rval = MBI->get_child_meshsets(newFace, children, 0); // all children
  MBERRCHK(rval, MBI);
  children.insert(newFace);
  rval = MBI->add_entities(newRootModel, children);
  MBERRCHK(rval, MBI);
  // add to the mentSelection of first model entity
  ModelEnt * firstModelEnt = (mentSelection.begin())->first;
  mentSelection[firstModelEnt].insert(newRootModel);
  return;
}

// construct the mesh: nothing to do, there is no mesh, really, only geometry
// in the form of mesh-based geometry
void MBSplitOp::execute_this()
{

}
}
