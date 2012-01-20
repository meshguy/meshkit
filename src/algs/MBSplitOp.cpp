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

// #include "meshkit/MBGeomOp.hpp"

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
  _globalId = -1;
  _direction[0]=_direction[1]=0.;
  _direction[2]=1.0;
  _closed = 1;
  _min_dot = 0.8;// acos(0.8)~= 36 degrees
}

MBSplitOp::~MBSplitOp()
{
}

//set up the crop/split of a mesh based surface
void MBSplitOp::set_options(int globalId, double dirx, double diry,
    double dirz, int closed, double min_dot)
{
  _globalId = globalId;
  /*for (int i=0; i<nPoints*3; i++)
   _polyline.push_back(polyline[i]);*/

  _direction[0] = dirx;
  _direction[1] = diry;
  _direction[2] = dirz;
  _closed = closed;
  _min_dot = min_dot;
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

}

// it is not a true mesh operation, it is a geometry operation
void MBSplitOp::execute_this()
{
  if (mentSelection.empty())
  {
    // change of strategy: do not return, go to
    // previous operation, find the result  (either another split or
    // an initial MBGeomOp, get the Range of the first mapped entity)
    // this will be the input to create other "ModelEnt"s, (mesh based geo)
    // find previous operation:
    GraphNode * prevNode = other_node(in_arcs());

    MeshOp * prevOp = reinterpret_cast<MeshOp*> (prevNode);
    std::string nameOp = prevNode->get_name();// just for debugging

    if (NULL==prevOp)
      return;

    MEntSelection & mentsel = prevOp->me_selection();
    // ments[0]
    ModelEnt * firstMe = (*(mentsel.begin()) ).first;
    moab::Range prevModelSet = mentsel[firstMe];

    if (prevModelSet.size()!=1)
      return;

    // this range has one set that can serve as basis for
    // mesh based geometry
    // the model ents have no model tag on the sets!!
    MEntVector model_ents; // existing model ents of dimension 2
    // the new ones will be put in the mentSelection!!!
    mk_core()->get_entities_by_dimension(2, model_ents);
    // assume moab_instance 0
    mk_core()->create_mbg_model_entities(prevModelSet[0], false);

    MEntVector model_ent_new; // existing model ents of dimension 2
        // the new ones will be put in the mentSelection!!!
    MEntVector model_ents_new; // all model ents of dimension 2...
    mk_core()->get_entities_by_dimension(2, model_ents_new);

    for (unsigned int i = model_ents.size(); i<model_ents_new.size(); i++)
    {
      moab::Range rr = mentSelection[model_ents_new[i]];
      rr.clear(); // just to use it , to avoid some warnings
    }
    //return;
  }

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
    MBERRCHK(rval, MBI);
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
  rval = pFacet->split_surface_with_direction(mset, _polyline, _direction, _closed, _min_dot,
      newFace);

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
}
