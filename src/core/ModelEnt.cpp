#include "ModelEnt.hpp"

namespace MeshKit
{
    
void ModelEnt::create_mesh_set(int flag) throw(Error)
{
  if (moabEntSet) return;

  if (!iGeomEnt) throw Error(MK_FAILURE, "Tried to create ModelEnt missing a geometry entity.");
  
    // need dimension, to tell whether to create a list or set
  if (-1 > flag || 1 < flag) throw Error(MK_FAILURE, "Invalid value for ordered flag.");
  
  moab::EntitySetProperty ordered = MESHSET_SET;
    // need type for later
  int this_tp = mkCore->iGeom.getEntType(iGeomEnt);
  if (1 == flag || (-1 == flag && iBase_EDGE == this_tp)) 
    ordered = MESHSET_ORDERED;
  
    // create the set
  moab::ErrorCode rval = mkCore->moabInstance->create_meshset(ordered, moabEntSet);
  MBCHKERR(rval, "Failed to create mesh set.");
  
    // tag it with this ModelEnt
  rval = mkCore->moabInstance->tag_set_data(mkCore->moabModelTag, &moabEntSet, 1, &this);
  MBERRCHK(rval, "Failed to set iMesh ModelEnt tag.");

    // relate the mesh to the geom
  err = mkCore->iRel.setEntSetRelation(mkCore->irel_pair(), iGeomEnt, iMeshEntSet);
  MKERRCHK(err, "Failed to set iRel relation for a mesh set.");

    // get parents and children, and link to corresponding sets; don't do this for any missing mesh sets, 
    // will get done when those sets get created
  std::vector<iGeom::EntityHandle> geom_ents;
  std::vector<iGeom::EntityHandle>::iterator vit;
  if (this_tp < iBase_REGION) {
    err = mkCore->iGeomInstance->getEntAdj(iGeomEnt, this_tp+1, geom_ents);
    MKERRCHK(err, "Trouble finding parent entities.");
    moab::EntityHandle seth;
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moabInstance->add_parent_child(moabEntSet, seth);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
  
  if (this_tp > iBase_VERTEX) {
    geom_ents.clear();
    err = mkCore->iGeomInstance->getEntAdj(iGeomEnt, this_tp-1, geom_ents);
    MKERRCHK(err, "Trouble finding child entities.");
    moab::EntityHandle seth;
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moabInstance->add_parent_child(seth, moabEntSet);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
}
    
int ModelEnt::dimension() const throw(Error) 
{
  iGeom_EntityType tp;
  iGeom::Error err = mkCore->iGeomInstance->getEntType(iGeomEnt, tp);
  IGERRCHK(err, "Couldn't get geom entity type.");
  return (tp - iBase_VERTEX);
}
  
double ModelEnt::measure() 
{
  double meas;
  iGeom::Error err = mkCore->iGeomInstance->measure(&iGeomEnt, 1, &meas);
  IGERRCHK(err, "Couldn't get geom entity measure.");
  return meas;
}
    
double double ModelEnt::measure_discrete() throw(Error)
{
  return measure();
}
      

  //- closest point to input
    double[3] ModelEnt::closest(double x, double y, double z) throw(Error)
{
  double closest[3];
  iGeom::Error err = mkCore->iGeomInstance->getEntClosestPt(iGeomEnt, x, y, z,
                                                            closest[0], closest[1], closest[2]);
  IGERRCHK(err, "Couldn't get geom entity closest point.");
  return closest;
}

  //- similar to closest_point, but based on mesh
double[3] ModelEnt::closest_discrete(double x, double y, double z) throw(Error) 
{
  return closest(x, y, z);
}
    

  //- return mesh bounding this model entity
double ModelEnt::boundary(int dim, std::vector<iMesh_EntityHandle> &mesh_ents);

}
