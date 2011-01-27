#include "iGeom.hh"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/Types.h"
#include "meshkit/VecUtil.hpp"
#include "moab/GeomTopoTool.hpp"

namespace MeshKit
{

ModelEnt::ModelEnt(MKCore *mk,
                   iGeom::EntityHandle geom_ent,
                   moab::EntityHandle mesh_ent,
                   int sizing_index) 
        : mkCore(mk), iGeomEnt(geom_ent), moabEntSet(mesh_ent), sizingFunctionIndex(sizing_index),
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
  std::vector<int> senses;
  moab::EntityHandle seth;
  moab::GeomTopoTool gt(mkCore->moab_instance(), false);
  if (this_tp < iBase_REGION) {
    err = mkCore->igeom_instance()->getEntAdj(iGeomEnt, (iBase_EntityType)(this_tp+1), geom_ents);
    IBERRCHK(err, "Trouble finding parent entities.");
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moab_instance()->add_parent_child(seth, moabEntSet);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
  
  if (this_tp > iBase_VERTEX) {
    geom_ents.clear();
    err = mkCore->igeom_instance()->getEntAdj(iGeomEnt, (iBase_EntityType)(this_tp-1), geom_ents);
    IBERRCHK(err, "Trouble finding child entities.");
    for (vit = geom_ents.begin(); vit != geom_ents.end(); vit++) {
      seth = mesh_handle(*vit);
      if (seth) {
        rval = mkCore->moab_instance()->add_parent_child(moabEntSet, seth);
        MBERRCHK(rval, "Trouble adding parent/child relation.");
      }
    }
  }
}

void ModelEnt::set_senses() 
{
  MEntVector adjs;
  int dim = dimension();
  moab::ErrorCode rval;
  iGeom::Error err;
  moab::GeomTopoTool gt(mkCore->moab_instance(), false);
  
    // downward, only for faces/regions
  if (1 < dim) {
    get_adjacencies(dim-1, adjs);
    
    for (MEntVector::iterator vit = adjs.begin(); vit != adjs.end(); vit++) {
      int dum_sense;
      if (2 == dim) (*vit)->set_upward_senses();
      else {
        err = mkCore->igeom_instance()->getEntNrmlSense((*vit)->geom_handle(), iGeomEnt, dum_sense);
        IBERRCHK(err, "Problem getting senses.");
        rval = gt.set_sense((*vit)->mesh_handle(), moabEntSet, dum_sense);
        MBERRCHK(rval, "Problem setting sense.");
      }
    }
  }

    // upward, only for edges/faces
  if (0 < dim && 3 > dim) {
    if (1 == dim) {
      set_upward_senses();
      return;
    }
    
    adjs.clear();
    get_adjacencies(dim+1, adjs);
    std::vector<moab::EntityHandle> ments;
    std::vector<int> senses;
    
    for (MEntVector::iterator vit = adjs.begin(); vit != adjs.end(); vit++) {
      int dum_sense;
      err = mkCore->igeom_instance()->getEntNrmlSense(iGeomEnt, (*vit)->geom_handle(), dum_sense);
      IBERRCHK(err, "Problem getting senses.");
      rval = gt.set_sense(moabEntSet, (*vit)->mesh_handle(), dum_sense);
      MBERRCHK(rval, "Problem setting sense.");
    }
  }
}

void ModelEnt::set_upward_senses() 
{
  assert("Should only call this for dimension == 1" && 1 == dimension());

  std::vector<ModelEnt*> adjs;
  std::vector<moab::EntityHandle> ments;
  std::vector<int> senses;
  get_adjacencies(2, adjs);
  int dum_sense;
  for (std::vector<ModelEnt*>::iterator vit = adjs.begin(); vit != adjs.end(); vit++) {
    iGeom::Error err = mkCore->igeom_instance()->getEgFcSense(iGeomEnt, (*vit)->geom_handle(), dum_sense);
    IBERRCHK(err, "Problem getting senses.");
    ments.push_back((*vit)->mesh_handle());
    senses.push_back(dum_sense);
  }
  moab::GeomTopoTool gt(mkCore->moab_instance());
  moab::ErrorCode rval = gt.set_senses(mesh_handle(), ments, senses);
  MBERRCHK(rval, "Problem setting senses on edge.");
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
  
void ModelEnt::get_adjacencies(int dim, MEntVector &adjs) const 
{
  std::vector<iGeom::EntityHandle> gents;
  iGeom::Error err = mkCore->igeom_instance()->getEntAdj(geom_handle(), (iBase_EntityType)dim, gents);
  IBERRCHK(err, "Trouble getting geom adjacencies.");
  adjs.resize(gents.size());
  err = mkCore->igeom_instance()->getArrData(&gents[0], adjs.size(), mkCore->igeom_model_tag(), 
                                             &adjs[0]);
  IBERRCHK(err, "Trouble getting ModelEnts for geom adjacencies.");
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

void ModelEnt::boundary(int dim,
                        MEntVector &entities,
                        std::vector<int> *senses,
                        std::vector<int> *group_sizes)
{
    // for a given surface, return the bounding edges in the form of edge groups,
    // oriented ccw around surface
  int result, ent_type;
  int this_dim;

    // shouldn't be calling this function if we're a vertex
  if (this_dim == iBase_VERTEX) throw Error(MK_BAD_INPUT, "Shouldn't call this for a vertex.");;
  
    // get adj ents & senses
  MEntVector tmp_ents;
  get_adjacencies(dim, tmp_ents);
  
    // special handling for edges and entities with single bounding entity - 
    // will just have 1 or 2 adj entities, set those directly and return
  if (1 == dim || 1 == tmp_ents.size()) {
    for (unsigned int i = 0; i < tmp_ents.size(); i++) {
      entities.push_back(tmp_ents[i]);
      if (group_sizes) group_sizes->push_back(1);
    }
    return;
  }
  
  MEntVector intersected_ents(tmp_ents.size());

    // get adjacent entities into a sorted, mutable list
  MEntVector b_ents;
  MEntSet dbl_curves, sgl_curves;
  b_ents = tmp_ents;
  std::sort(b_ents.begin(), b_ents.end());

  while (!b_ents.empty()) {
      // get 1st in group
    
    MEntVector group_stack, this_group;
    group_stack.push_back(b_ents.front());
    
      // while there are still entities on the stack
    while (!group_stack.empty()) {
      
        // pop one off & get its sense
      ModelEnt *this_entity = group_stack.back(); group_stack.pop_back();
      
      int this_sense;
      iGeom::Error err = 
          mkCore->igeom_instance()->getSense(this_entity->geom_handle(), geom_handle(), this_sense);

        // if we already have this one, continue
      if ((0 != this_sense && 
           std::find(this_group.begin(), this_group.end(), this_entity) != this_group.end()) ||
          std::find(dbl_curves.begin(), dbl_curves.end(), this_entity) != dbl_curves.end())
        continue;

        // either way we need the d-2 entities
      MEntVector bridges;
      get_adjacencies(dim-2, bridges);
      
        // only remove from the list of candidates if it's not dual-sensed
      if (0 != this_sense) 
        b_ents.erase(std::remove(b_ents.begin(), b_ents.end(), this_entity), b_ents.end());

        // if it's double-sensed and this is the first time we're seeing it, find the right sense
      if (0 == this_sense) {
        if (std::find(sgl_curves.begin(), sgl_curves.end(), this_entity) == sgl_curves.end()) {
          sgl_curves.insert(this_entity);
          if (!this_group.empty()) {
            ModelEnt *common_v = this_entity->shared_entity(this_group.back(), 0);
            if (common_v == bridges[0]) this_sense = 1;
            else if (common_v == bridges[1]) this_sense = -1;
            else return;
          }
            // else, if this is the first one in the loop, just choose a sense
          else {
            this_sense = 1;
          }
        }
        else {
            // else this is the second time we're seeing it, move it to the dbl_curves list and remove
            // from the candidates list
          dbl_curves.insert(this_entity);
          sgl_curves.erase(this_entity);
          b_ents.erase(std::remove(b_ents.begin(), b_ents.end(), this_entity), b_ents.end());
        }
      }
          
        // it's in the group; put on group & remove from untreated ones
      this_group.push_back(this_entity);
      if (senses) senses->push_back(this_sense);
      MEntVector tmp_from, tmp_adjs;

        // if we're on a face and we're the first in a group, check sense of this first
        // edge; make sure "next" in loop sense is last on list
      if (2 == dim && this_group.size() == 1) {

          // get vertex which we know is shared by the "right" next edge; first get the vertices
          // if sense of current edge is forward, it's the 2nd vertex we want,
          // otherwise the first
        if (1 == this_sense && bridges.size() > 1) tmp_from.push_back(bridges[1]);
        else tmp_from.push_back(bridges[0]);
        tmp_adjs = b_ents;
        get_adjs_bool(tmp_from, dim-1, tmp_adjs, INTERSECT);
      }
      else {
        std::copy(bridges.begin(), bridges.end(), std::back_inserter(tmp_from));
        MEntVector tmp_adjs2;
        get_adjs_bool(tmp_from, dim-1, tmp_adjs2, UNION);
        std::sort(tmp_adjs2.begin(), tmp_adjs2.end());
        tmp_adjs.resize(tmp_adjs2.size());
        tmp_adjs.erase(std::set_intersection(b_ents.begin(), b_ents.end(),
                                             tmp_adjs2.begin(), tmp_adjs2.end(),
                                             tmp_adjs.begin()), tmp_adjs.end());
      }

      if (2 == dim && tmp_adjs.size() > 1) {
          // more than one adjacent edge - need to evaluate winding angle to find right one
        ModelEnt *next_ent = next_winding(this_entity, this_sense, tmp_adjs);
        if (NULL == next_ent) return;
        group_stack.push_back(next_ent);
      }
      else if (!tmp_adjs.empty()) 
        std::copy(tmp_adjs.begin(), tmp_adjs.end(), std::back_inserter(group_stack));
    }
    
    // put group in group list
    std::copy(this_group.begin(), this_group.end(), std::back_inserter(entities));
    if (group_sizes) group_sizes->push_back(this_group.size());
  }
}

ModelEnt *ModelEnt::shared_entity(ModelEnt *ent2, int to_dim) 
{
    // find the shared entity between the two entities of the prescribed dimension
  MEntVector from_ents(2), to_ents;
  from_ents[0] = this;
  from_ents[1] = ent2;
  try {
    get_adjs_bool(from_ents, to_dim, to_ents, INTERSECT); 
  }
  catch (Error err) {
    if (err.error_code() != MK_SUCCESS || to_ents.empty()) return NULL;
  }

  if (to_ents.size() > 1) throw Error(MK_MULTIPLE_FOUND, "Multiple shared entities found.");
  
  return to_ents[0];
}

void ModelEnt::get_adjs_bool(MEntVector &from_ents,
                             int to_dim,
                             MEntVector &to_ents,
                             BooleanType op_type) 
{
  if (from_ents.empty()) {
    to_ents.clear();
    return;
  }

  int result;
  MEntVector bridges;

  MEntVector::iterator from_it = from_ents.begin();
  if (to_ents.empty() && op_type == INTERSECT) {
    from_ents.front()->get_adjacencies(to_dim, bridges);
    std::copy(bridges.begin(), bridges.end(), std::back_inserter(to_ents));
    from_it++;
  }

  std::sort(to_ents.begin(), to_ents.end());
  MEntVector result_ents(to_ents.size());
  for (; from_it != from_ents.end(); from_it++) {
    bridges.clear();
    (*from_it)->get_adjacencies(to_dim, bridges);
    if (op_type == INTERSECT) {
      std::sort(bridges.begin(), bridges.end());
      result_ents.erase(std::set_intersection(to_ents.begin(), to_ents.end(),
                                              bridges.begin(), bridges.end(),
                                              result_ents.begin()), result_ents.end());
    }
    else {
      std::copy(bridges.begin(), bridges.end(), std::back_inserter(result_ents));
    }
    
    to_ents = result_ents;
  }
  
  if (op_type == UNION)
    std::sort(to_ents.begin(), to_ents.end());
  
  to_ents.erase(std::unique(to_ents.begin(), to_ents.end()), to_ents.end());
}

ModelEnt *ModelEnt::next_winding(ModelEnt *this_edge, 
                                 int this_sense, 
                                 MEntVector &tmp_adjs) 
{
    // given this_entity, a d-1 entity bounding gentity, and optional "next"
    // entities in tmp_adjs, find the next one based on windings around the shared
    // vertex
  
    // first, get the shared vertex 
  MEntVector verts;
  this_edge->get_adjacencies(0, verts);

  ModelEnt *shared_vert = verts[0];
  if (this_sense == 1 && verts.size() > 1) shared_vert = verts[1];

    // get locations just before the vertex, at the vertex, and just after the vertex
  double v1[3], v2[3], v3[3];
  double umin, umax;
  mkCore->igeom_instance()->getEntURange(this_edge->geom_handle(), umin, umax);
  double utgt;
  if (1 == this_sense) utgt = umin + 0.9 * (umax - umin);
  else utgt = umin + 0.1 * (umax - umin);
  mkCore->igeom_instance()->getEntUtoXYZ(this_edge->geom_handle(), utgt, v1[0], v1[1], v1[2]);
  mkCore->igeom_instance()->getVtxCoord(shared_vert->geom_handle(), v2[0], v2[1], v2[2]);
  double v21[] = {v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]};
  VecUtil::normalize(v21);

    // get the normal vector at the vertex
  double normal[3];
  mkCore->igeom_instance()->getEntNrmlXYZ(geom_handle(), v2[0], v2[1], v2[2],
                                          normal[0], normal[1], normal[2]);
  
    // now loop over candidates, finding magnitude of swept angle
  ModelEnt *other = NULL;
  double angle = VecUtil::TWO_PI;
  bool is_adj;
  for (MEntVector::iterator vit = tmp_adjs.begin(); vit != tmp_adjs.end(); vit++) {
      // if we're here, we have multiple candidates, therefore don't choose the same one
    if (*vit == this_edge)
      continue;
    mkCore->igeom_instance()->isEntAdj((*vit)->geom_handle(), shared_vert->geom_handle(), is_adj);
    if (!is_adj)
      continue;

      // get param range
    mkCore->igeom_instance()->getEntURange((*vit)->geom_handle(), umin, umax);
      // get sense
    int tmp_sense;
    mkCore->igeom_instance()->getSense((*vit)->geom_handle(), geom_handle(), tmp_sense);
    if (1 == tmp_sense) utgt = umin + 0.1 * (umax - umin);
    else utgt = umin + 0.9 * (umax - umin);
    mkCore->igeom_instance()->getEntUtoXYZ((*vit)->geom_handle(), utgt, v3[0], v3[1], v3[2]);
    double v23[] = {v3[0]-v2[0], v3[1]-v2[1], v3[2]-v2[2]};
    VecUtil::normalize(v23);

    VecUtil::cross(normal, v23, v1);
    
    VecUtil::cross(v1, normal, v3);

    double x = VecUtil::dot(v21, v3);
    double y = VecUtil::dot(v21, v1);

    assert(x != 0.0 || y != 0.0);
    double this_angle = atan2(y, x);
    if (this_angle < 0.0)
    {
      this_angle += VecUtil::TWO_PI;
    }
    if (this_angle < angle) {
      other = *vit;
      angle = this_angle;
    }
  }
  
  return other;
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

moab::EntityHandle ModelEnt::mesh_handle(iGeom::EntityHandle gent) const 
{
    // use iRel to get this information
  moab::EntityHandle ment = 0;
  Error err = mkCore->irel_pair()->getEntSetRelation(gent, false, IBSHR(ment));
  MKERRCHK(err, "Failed to get mesh set handle for geometry entity.");
  return ment;
}

} // namespace MeshKit
