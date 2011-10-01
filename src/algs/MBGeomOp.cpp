/*
 * MBGeomOp.cpp
 *
 *  Created on: Sep 30, 2011
 */

#include <iostream>

#include "meshkit/MBGeomOp.hpp"
#include "meshkit/ModelEnt.hpp"

#include "moab/Core.hpp"
#include "moab/GeomTopoTool.hpp"


namespace MeshKit {

//Entity Type initialization for geometrization; no mesh output
moab::EntityType MBGeomOp_tps[] = { moab::MBMAXTYPE }; // no mesh, really
const moab::EntityType* MBGeomOp::output_types()
{
  return MBGeomOp_tps;
}

MBGeomOp::MBGeomOp(MKCore *mk_core, const MEntVector &me_vec) :
      MeshScheme(mk_core, me_vec)
{
  // TODO Auto-generated constructor stub

}

MBGeomOp::~MBGeomOp()
{
  // TODO Auto-generated destructor stub
}

//set up the geometrization of a model ent (face) for mesh-based geometry
void MBGeomOp::setup_this()
{
  if (mentSelection.empty())
    return;
  ModelEnt * me = mentSelection.begin()->first;
  // get the one and only model entity, which should be of dimension 2, and geometrize it
  // this should be the model entity handle
  moab::EntityHandle mset = me->mesh_handle();
  iGeom::EntityHandle gent = me->geom_handle() ;
  if (gent!=0)
  {
    std::cout<< "geometry exists on the model ent, do not geometrize\n";
    return; // we do have something here, do not "geometrize"
  }
  moab::GeomTopoTool gtt(mk_core()->moab_instance());
  moab::EntityHandle output;
  moab::ErrorCode rval = gtt.geometrize_surface_set(mset, output);
  MBERRCHK(rval, mk_core()->moab_instance());
  // now, get all entity sets in here; the one with dimension 2, put it in mset, and
  // form other model ents with other sets;
  moab::Range geom_sets;
  moab::Tag geomTag=mk_core()->moab_geom_dim_tag();

  rval = mk_core()->moab_instance()->get_entities_by_type_and_tag(output, moab::MBENTITYSET, &geomTag, NULL, 1, geom_sets,
            moab::Interface::UNION);
  MBERRCHK(rval, mk_core()->moab_instance());

  std::vector <int> dims;
  dims.resize(geom_sets.size());
  rval = mk_core()->moab_instance()->tag_get_data(geomTag, geom_sets, &dims[0]);
  MBERRCHK(rval, mk_core()->moab_instance());

  // find the one of dimension 2, it is the main one:
  // the rest should be of dimension 1 or 0, and for each we should have a ModelEnt
  // we should maybe repeat the process from populate mesh...
  // maybe we should create new model ents for these new gents; they will be meshable with
  //  camal, if the geo handles are from FBiGeom

  return;

  // get dimensions

}

  // construct the mesh: nothing to do, there is no mesh, really, only geometry
  // in the form of mesh-based geometry
void MBGeomOp::execute_this()
{

}

}
