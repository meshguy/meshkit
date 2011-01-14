#include "meshkit/ModelEnt.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/Types.h"

namespace MeshKit
{

ModelEnt::ModelEnt(MKCore *mk,
                   iGeom::EntityHandle geom_ent,
                   moab::EntityHandle mesh_ent,
                   const SizingFunction &mesh_size) 
        : mkCore(mk), iGeomEnt(geom_ent), moabEntSet(mesh_ent), sizingFunction(mesh_size),
          meshIntervals(-1), intervalFirmness(SOFT), meshedState(NO_MESH) 
{}

ModelEnt::~ModelEnt() 
{}

void ModelEnt::create_mesh_set(int flag)
{
  if (moabEntSet) return;

  if (!iGeomEnt) throw Error(MK_FAILURE, "Tried to create ModelEnt missing a geometry entity.");
  
    // need dimension, to tell whether to create a list or set
  if (-1 > flag || 1 < flag) throw Error(MK_FAILURE, "Invalid value for ordered flag.");
  
  moab::EntitySetProperty ordered = moab::MESHSET_SET;
    // need type for later
  iBase_EntityType this_tp;
  iGeom::Error err = mkCore->igeom_instance()->getEntType(iGeomEnt, this_tp);
  IBERRCHK(err, "Trouble getting entity type.");
  if (1 == flag || (-1 == flag && iBase_EDGE == this_tp)) 
    ordered = moab::MESHSET_ORDERED;
  
    // create the set
  moab::ErrorCode rval = mkCore->moab_instance()->create_meshset(ordered, moabEntSet);
  MBERRCHK(rval, "Failed to create mesh set.");
  
    // tag it with this ModelEnt
  ModelEnt *this_ptr = this;
  rval = mkCore->moab_instance()->tag_set_data(mkCore->moab_model_tag(), &moabEntSet, 1, &this_ptr);
  MBERRCHK(rval, "Failed to set iMesh ModelEnt tag.");

    // relate the mesh to the geom
  Error merr = mkCore->irel_pair()->setEntSetRelation(iGeomEnt, IBSH(moabEntSet));
  MKERRCHK(merr, "Failed to set iRel relation for a mesh set.");

    // get parents and children, and link to corresponding sets; don't do this for any missing mesh sets, 
    // will get done when those sets get created
  std::vector<iGeom::EntityHandle> geom_ents;
  std::vector<iGeom::EntityHandle>::iterator vit;
  if (this_tp < iBase_REGION) {
    err = mkCore->igeom_instance()->getEntAdj(iGeomEnt, (iBase_EntityType)(this_tp+1), geom_ents);
    IBERRCHK(err, "Trouble finding parent entities.");
    moab::EntityHandle seth;
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moab_instance()->add_parent_child(moabEntSet, seth);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
  
  if (this_tp > iBase_VERTEX) {
    geom_ents.clear();
    err = mkCore->igeom_instance()->getEntAdj(iGeomEnt, (iBase_EntityType)(this_tp-1), geom_ents);
    IBERRCHK(err, "Trouble finding child entities.");
    moab::EntityHandle seth;
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moab_instance()->add_parent_child(seth, moabEntSet);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
}
    
/** \brief Commit mesh to a model entity
 *
 * Takes the input mesh entities, adds them to the entity set for this model entity,
 * and (if both-type relation on the mesh side) sets the relations to the corresponding
 * geometry entity.
 * \param mesh_ents Mesh entities being assigned to this model entity
 * \param mstate The meshed state after this mesh is added
 */
void ModelEnt::commit_mesh(moab::Range &mesh_ents, MeshedState mstate) 
{
  std::vector<moab::EntityHandle> ent_vec;
  std::copy(mesh_ents.begin(), mesh_ents.end(), ent_vec.begin());
  commit_mesh(ent_vec, mstate);
}

/** \brief Commit mesh to a model entity
 *
 * Takes the input mesh entities, adds them to the entity set for this model entity,
 * and (if both-type relation on the mesh side) sets the relations to the corresponding
 * geometry entity.
 * \param mesh_ents Mesh entities being assigned to this model entity
 * \param mstate The meshed state after this mesh is added
 */
void ModelEnt::commit_mesh(std::vector<moab::EntityHandle> &mesh_ents,
                           MeshedState mstate) 
{
  iRel::Error err =
      mkCore->irel_pair()->setEntEntArrRelation(geom_handle(), false,
                                                (iBase_EntityHandle*)&mesh_ents[0], mesh_ents.size());
  IBERRCHK(err, "Trouble committing mesh.");

  meshedState = mstate;
}
  
int ModelEnt::dimension() const 
{
  iBase_EntityType tp;
  iGeom::Error err = mkCore->igeom_instance()->getEntType(iGeomEnt, tp);
  IBERRCHK(err, "Couldn't get geom entity type.");
  return (tp - iBase_VERTEX);
}
  
double ModelEnt::measure() const
{
  double meas;
  iGeom::Error err = mkCore->igeom_instance()->measure(&iGeomEnt, 1, &meas);
  IBERRCHK(err, "Couldn't get geom entity measure.");
  return meas;
}
    
double ModelEnt::measure_discrete() const
{
  return measure();
}

  //- closest point to input
void ModelEnt::closest(double x, double y, double z, double *close) const
{
  iGeom::Error err = mkCore->igeom_instance()->getEntClosestPt(iGeomEnt, x, y, z,
                                                            close[0], close[1], close[2]);
  IBERRCHK(err, "Couldn't get geom entity closest point.");
}

  //- similar to closest_point, but based on mesh
void ModelEnt::closest_discrete(double x, double y, double z, double *close) const
{
  return closest(x, y, z, close);
}
    
  //- return mesh bounding this model entity
void ModelEnt::boundary(int dim, std::vector<moab::EntityHandle> &mesh_ents) const
{
  throw Error(MK_NOT_IMPLEMENTED, "Not implemented yet.");
}

void ModelEnt::boundary(int dim, 
                        std::vector<moab::EntityHandle> &forward_ents,
                        std::vector<moab::EntityHandle> &reverse_ents) const 
{
  throw Error(MK_NOT_IMPLEMENTED, "Not implemented yet.");
}

void ModelEnt::boundary(int dim, 
                        moab::Range &ents) const
{
  throw Error(MK_NOT_IMPLEMENTED, "Not implemented yet.");
}

iGeom::EntityHandle ModelEnt::geom_handle(moab::EntityHandle ment) const 
{
    // use iRel to get this information
  iBase_EntityHandle gent = 0;
  Error err = mkCore->irel_pair()->getSetEntRelation(IBSH(ment), true, gent);
  MKERRCHK(err, "Failed to get geometry handle for mesh set.");
  return gent;
}

moab::EntityHandle ModelEnt::mesh_handle(iBase_EntityHandle gent) const 
{
    // use iRel to get this information
  moab::EntityHandle ment = 0;
  Error err = mkCore->irel_pair()->getEntSetRelation(gent, false, IBSHR(ment));
  MKERRCHK(err, "Failed to get mesh set handle for geometry entity.");
  return ment;
}

} // namespace MeshKit
