#include "meshkit/iGeom.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/Types.hpp"
#include "meshkit/VecUtil.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/CN.hpp"
#include "RefEntity.hpp"

namespace MeshKit
{

ModelEnt::ModelEnt(MKCore *mk,
                   iGeom::EntityHandle geom_ent,
                   moab::EntityHandle mesh_ent,
                   int sizing_index) 
        : mkCore(mk), iGeomEnt(geom_ent), iGeomSet(NULL), moabEntSet(mesh_ent), sizingFunctionIndex(sizing_index),
          meshIntervals(-1), intervalFirmness(DEFAULT), meshedState(NO_MESH) 
{}

ModelEnt::ModelEnt(MKCore *mk,
                   iGeom::EntitySetHandle geom_ent,
                   moab::EntityHandle mesh_ent,
                   int sizing_index) 
        : mkCore(mk), iGeomEnt(NULL), iGeomSet(geom_ent), moabEntSet(mesh_ent), sizingFunctionIndex(sizing_index),
          meshIntervals(-1), intervalFirmness(DEFAULT), meshedState(NO_MESH) 
{}

ModelEnt::~ModelEnt() 
{}

void ModelEnt::sizing_function_index(int ind, bool children_too)
{
  sizingFunctionIndex = ind;

    // set on children too, if requested
  if (children_too && dimension() > 1) {
    MEntVector childrn;
    get_adjacencies(dimension()-1, childrn);
    for (MEntVector::iterator vit = childrn.begin(); vit != childrn.end(); vit++) 
      (*vit)->sizing_function_index(ind, children_too);
  }
}
    
void ModelEnt::create_mesh_set(int flag)
{
  if (moabEntSet) return;

  if (!iGeomEnt && !iGeomSet) 
    throw Error(MK_FAILURE, "Tried to create ModelEnt missing a geometry entity or set.");
  
    // need dimension, to tell whether to create a list or set
  if (-1 > flag || 1 < flag) throw Error(MK_FAILURE, "Invalid value for ordered flag.");
  
  moab::EntitySetProperty ordered = moab::MESHSET_SET;
    // need type for later
  iBase_EntityType this_tp;
  iGeom::Error err;
  if (iGeomEnt) {
    err = mkCore->igeom_instance()->getEntType(iGeomEnt, this_tp);
    IBERRCHK(err, "Trouble getting entity type.");
    if (1 == flag || (-1 == flag && iBase_EDGE == this_tp)) 
      ordered = moab::MESHSET_ORDERED;
  }
  
    // create the set
  moab::ErrorCode rval = mkCore->moab_instance()->create_meshset(ordered, moabEntSet);
  MBERRCHK(rval, mkCore->moab_instance());
  
    // tag it with this ModelEnt
  ModelEnt *this_ptr = this;
  rval = mkCore->moab_instance()->tag_set_data(mkCore->moab_model_tag(), &moabEntSet, 1, &this_ptr);
  MBERRCHK(rval, mkCore->moab_instance());

    // relate the mesh to the geom
  iRel::Error merr;
  if (iGeomEnt) {
    merr = mkCore->irel_pair()->setEntSetRelation(iGeomEnt, IBSH(moabEntSet));
  }
  else {
    merr = mkCore->group_set_pair()->setSetSetRelation(iGeomSet, IBSH(moabEntSet));
  }
  IBERRCHK(merr, "Failed to set iRel relation for a mesh set.");

    // if this is a group, we're done
  if (iGeomSet) return;

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
      try {seth = mesh_handle(*vit);
      }
      catch (Error err) {
            // just continue, will get this when parent gets created later
          continue;
        }
      if (seth) {
        rval = mkCore->moab_instance()->add_parent_child(seth, moabEntSet);
        MBERRCHK(rval, mkCore->moab_instance());
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
        MBERRCHK(rval, mkCore->moab_instance());
      }
    }
  }
}

// find any existing set to relate
bool ModelEnt::exist_mesh_set()
{
  if (!iGeomEnt) return false;

  int dim = dimension();
  int id = reinterpret_cast<RefEntity*> (iGeomEnt)->id();
  const void* dims[] = {&dim};

  moab::Range ents;
  moab::Tag geom_dim_tag = mkCore->moab_geom_dim_tag();
  moab::ErrorCode rval = mkCore->moab_instance()->
    get_entities_by_type_and_tag(0, moab::MBENTITYSET, &geom_dim_tag,
				 dims, 1, ents);
  MBERRCHK(rval, mkCore->moab_instance());

  moab::EntityHandle meshset;
  rval = mkCore->moab_instance()->create_meshset(moab::MESHSET_ORDERED, meshset);
  MBERRCHK(rval, mkCore->moab_instance());

  rval = mkCore->moab_instance()->add_entities(meshset, ents);
  MBERRCHK(rval, mkCore->moab_instance());

  const void* ids[] = {&id};
  ents.clear();
  moab::Tag id_tag = mkCore->moab_global_id_tag();
  rval = mkCore->moab_instance()->
    get_entities_by_type_and_tag(meshset, moab::MBENTITYSET, &id_tag,
				 ids, 1, ents);
  MBERRCHK(rval, mkCore->moab_instance());

  if (!ents.size()) return false;

  moabEntSet = ents[0];
  iRel::Error merr = mkCore->irel_pair()->setEntSetRelation(iGeomEnt, IBSH(moabEntSet));
  IBERRCHK(merr, "Failed to set iRel relation for a mesh set.");

  return true;
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
        MBERRCHK(rval, mkCore->moab_instance());
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
      MBERRCHK(rval, mkCore->moab_instance());
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
  MBERRCHK(rval, mkCore->moab_instance());
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
  std::copy(mesh_ents.begin(), mesh_ents.end(), std::back_inserter(ent_vec));
  commit_mesh(&ent_vec[0], ent_vec.size(), mstate);
}

/** \brief Commit mesh to a model entity
 *
 * Takes the input mesh entities, adds them to the entity set for this model entity,
 * and (if both-type relation on the mesh side) sets the relations to the corresponding
 * geometry entity.
 * \param mesh_ents Mesh entities being assigned to this model entity
 * \param mstate The meshed state after this mesh is added
 */
void ModelEnt::commit_mesh(moab::EntityHandle *mesh_ents,
                           int num_ents,
                           MeshedState mstate) 
{
    // if it's BOTH, set the relation through iRel
  if (mkCore->irel_pair()->get_rel_type(1) == iRel::BOTH) {
    iRel::Error err =
        mkCore->irel_pair()->setEntEntArrRelation(geom_handle(), false,
                                                  (iBase_EntityHandle*)mesh_ents, num_ents);
    IBERRCHK(err, "Trouble committing mesh.");
  }

    // either way, add them to the set
  moab::ErrorCode rval = mkCore->moab_instance()->add_entities(moabEntSet, mesh_ents, num_ents);
  MBERRCHK(rval, mkCore->moab_instance());
  
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

void ModelEnt::evaluate(double x, double y, double z, 
                        double *close,
                        double *direction,
                        double *curvature1,
                        double *curvature2) const
{
  iGeom::Error err;
  if (0 == dimension() || 3 == dimension()) {
    if (direction || curvature1 || curvature2) {
      MKERRCHK(Error(MK_BAD_INPUT), "Direction or curvature not available for entities of this type.");
    }
    else if (!close) {
      MKERRCHK(Error(MK_BAD_INPUT), "Must pass closest point pointer for output.");
    }
  }

  if (0 == dimension()) {
    err = mk_core()->igeom_instance()->getVtxCoord(geom_handle(), close[0], close[1], close[2]);
    IBERRCHK(err, "Problem getting vertex coordinates.");
    return;
  }
  else if (3 == dimension()) {
    err = mk_core()->igeom_instance()->getEntClosestPt(geom_handle(), x, y, z, 
                                                       close[0], close[1], close[2]);
    IBERRCHK(err, "Problem getting closest point for this model entity.");
  }
  else {
    double cls[3], dir[3], curve1[3], curve2[3];
    if (!close) close = cls;
    if (!direction) direction = dir;
    if (!curvature1) curvature1 = curve1;
    if (1 == dimension()) 
      err = mk_core()->igeom_instance()->getEgEvalXYZ(geom_handle(), x, y, z,
                                                      close[0], close[1], close[2],
                                                      direction[0], direction[1], direction[2], 
                                                      curvature1[0], curvature1[1], curvature1[2]);
    else if (2 == dimension()) {
      if (!curvature2) curvature2 = curve2;
      err = mk_core()->igeom_instance()->getFcEvalXYZ(geom_handle(), x, y, z,
                                                      close[0], close[1], close[2],
                                                      direction[0], direction[1], direction[2], 
                                                      curvature1[0], curvature1[1], curvature1[2], 
                                                      curvature2[0], curvature2[1], curvature2[2]);
    }

    IBERRCHK(err, "Couldn't evaluate model entity.");
  }
}

int ModelEnt::id() const
{
    // get a global id for this entity, if there is one
  if (geom_handle()) {
    iGeom::TagHandle gid_tag;
    iGeom::Error err = mk_core()->igeom_instance()->getTagHandle("GLOBAL_ID", gid_tag);
    if (iBase_SUCCESS == err) {
      int gid;
      err = mk_core()->igeom_instance()->getIntData(geom_handle(), gid_tag, gid);
      if (iBase_SUCCESS == err) return gid;
    }
  }
  
  if (mesh_handle()) {
    moab::Tag gid_tag;
    moab::ErrorCode err = mk_core()->moab_instance()->tag_get_handle("GLOBAL_ID", sizeof(int), moab::MB_TAG_SPARSE, 
                                                                     moab::MB_TYPE_INTEGER, gid_tag);
    if (moab::MB_SUCCESS == err) {
      int gid;
      moab::EntityHandle this_mesh = mesh_handle();
      err = mk_core()->moab_instance()->tag_get_data(gid_tag, &this_mesh, 1, &gid);
      if (moab::MB_SUCCESS == err) return gid;
    }
  }
  
  return -1;
}

void ModelEnt::get_mesh(int dim,
                        std::vector<moab::EntityHandle> &ments,
                        bool bdy_too)
{
  if (dim > dimension()) throw Error(MK_BAD_INPUT, "Called get_mesh for dimension ", dim, " on a ", 
                                     dimension(), "-dimensional entity.");

  if (!bdy_too || dim == dimension()) {
      // just owned entities, which will be in the set
    moab::ErrorCode rval = mk_core()->moab_instance()->get_entities_by_dimension(moabEntSet, dimension(), ments);
    MBERRCHK(rval, mkCore->moab_instance());
    return;
  }

    // filter out the cases where order matters and boundary is requested too
  else if (1 == dimension() && 0 == dim) {
    std::vector<moab::EntityHandle> tmp_edgevs, tmp_vvs;
    moab::ErrorCode rval = mk_core()->moab_instance()->get_entities_by_dimension(moabEntSet, 0, tmp_edgevs);
    MBERRCHK(rval, mkCore->moab_instance());
    boundary(0, tmp_vvs);
    if (tmp_vvs.empty()) throw Error(MK_NOT_FOUND, "No mesh on bounding entity.");
    ments.push_back(tmp_vvs[0]);
    std::copy(tmp_edgevs.begin(), tmp_edgevs.end(), std::back_inserter(ments));
    if (2 == tmp_vvs.size()) ments.push_back(tmp_vvs[1]);
    else if (1 == tmp_vvs.size()) ments.push_back(tmp_vvs[0]);
    return;
  }
  else {
      // remaining case is bdy_too=true and dim < dimension()
      // first, get dimension()-dimensional entities, then get dim-dimensional neighbors of those
    std::vector<moab::EntityHandle> tmp_ments;
    moab::ErrorCode rval = mk_core()->moab_instance()->get_entities_by_dimension(moabEntSet, dimension(), tmp_ments);
    MBERRCHK(rval, mkCore->moab_instance());
    rval = mk_core()->moab_instance()->get_adjacencies(&tmp_ments[0], tmp_ments.size(), dim, false, ments, 
                                                       moab::Interface::UNION);
    MBERRCHK(rval, mkCore->moab_instance());
  }
}
    
void ModelEnt::boundary(int dim,
                        std::vector<moab::EntityHandle> &bdy,
                        std::vector<int> *senses,
                        std::vector<int> *group_sizes) 
{
  MEntVector me_loop;
  std::vector<int> me_senses, me_group_sizes;
  boundary(dim, me_loop, &me_senses, &me_group_sizes);
  
  MEntVector::iterator vit = me_loop.begin();
  std::vector<int>::iterator sense_it = me_senses.begin(), grpsize_it = me_group_sizes.begin();
  int gents_ctr = 0, first_ment = 0;

    // packing into a single list, but need to keep track of loop/shell sizes; don't increment
    // grpsize_it here, just when we've done all the mes in this group
  for (vit = me_loop.begin(); vit != me_loop.end(); vit++, sense_it++) {
    std::vector<moab::EntityHandle> tmp_ments;
    (*vit)->get_mesh(dim, tmp_ments, true);
    if (*sense_it == SENSE_REVERSE) 
      std::reverse(tmp_ments.begin(), tmp_ments.end());
    if (2 == dimension() && 0 == dim)
        // assembling vertices on loops; don't do last, since that'll be the first of the
        // next edge
      std::copy(tmp_ments.begin(), tmp_ments.end()-1, std::back_inserter(bdy));
    else
      std::copy(tmp_ments.begin(), tmp_ments.end(), std::back_inserter(bdy));
    if (senses) {
      for (unsigned int i = 0; i < tmp_ments.size(); i++) senses->push_back(*sense_it);
    }
    if (group_sizes) {
      gents_ctr++;
      if (gents_ctr == *grpsize_it) {
        group_sizes->push_back(bdy.size()-first_ment);
          // first_ment will be the first element in the next group
        first_ment = bdy.size();
        gents_ctr = 0;
        grpsize_it++;
      } // block at end of a given loop
    } // block updating group size
  } // over loop
  
}
    
void ModelEnt::boundary(int dim,
                        MEntVector &entities,
                        std::vector<int> *senses,
                        std::vector<int> *group_sizes)
{
    // for a given entity, return the bounding edges in the form of edge groups,
    // oriented ccw around surface
  int result, ent_type;
  int this_dim = dimension();

    // shouldn't be calling this function if we're a vertex
  if (this_dim == iBase_VERTEX) throw Error(MK_BAD_INPUT, "Shouldn't call this for a vertex.");
  else if (dim >= this_dim) throw Error(MK_BAD_INPUT, "Calling for boundary entities of equal or greater dimension.");;
  
    // get (d-1)-adj ents & senses
  MEntVector tmp_ents;
  get_adjacencies(this_dim-1, tmp_ents);

    // get out here if we're a curve
  if (1 == this_dim) {
    for (unsigned int i = 0; i < tmp_ents.size(); i++) {
      entities.push_back(tmp_ents[i]);
      if (senses) senses->push_back(SENSE_FORWARD);
    }
    return;
  }
  
    // get adjacent entities into a sorted, mutable list
  MEntVector b_ents = tmp_ents;
  MEntSet dbl_curves, sgl_curves;
  
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
      this_entity->get_adjacencies(dimension()-2, bridges);
      
        // only remove from the list of candidates if it's not dual-sensed
      if (0 != this_sense) 
        b_ents.erase(std::remove(b_ents.begin(), b_ents.end(), this_entity), b_ents.end());

        // if it's double-sensed and this is the first time we're seeing it, find the right sense
      else {
        assert(dimension() == 2);
        if (std::find(sgl_curves.begin(), sgl_curves.end(), this_entity) == sgl_curves.end()) {
          sgl_curves.insert(this_entity);
          if (!this_group.empty()) {
            ModelEnt *common_v = this_entity->shared_entity(this_group.back(), 0);
            if (common_v == bridges[0]) this_sense = 1;
            else if (common_v == bridges[1]) this_sense = -1;
            else assert(false);
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
          if (senses) {
              // need to get the sense of the first appearance in list
            MEntVector::iterator vit = std::find(this_group.begin(), this_group.end(), this_entity);
            assert(vit != this_group.end());
            int other_sense = (*senses)[vit-this_group.begin()];
            this_sense = SENSE_REVERSE*other_sense;
          }
        }
      }
          
        // it's in the group; put on group & remove from untreated ones
      this_group.push_back(this_entity);
      if (senses) senses->push_back(this_sense);
      MEntVector tmp_from, tmp_adjs;

        // if we're on a face and we're the first in a group, check sense of this first
        // edge; make sure "next" in loop sense is last on list
      if (2 == dimension() && this_group.size() == 1) {

          // get vertex which we know is shared by the "right" next edge; first get the vertices
          // if sense of current edge is forward, it's the 2nd vertex we want,
          // otherwise the first
        if (1 == this_sense && bridges.size() > 1) tmp_from.push_back(bridges[1]);
        else tmp_from.push_back(bridges[0]);
        tmp_adjs = b_ents;
        get_adjs_bool(tmp_from, dimension()-1, tmp_adjs, INTERSECT);
      }
      else {
        std::copy(bridges.begin(), bridges.end(), std::back_inserter(tmp_from));
        MEntVector tmp_adjs2;
        get_adjs_bool(tmp_from, dimension()-1, tmp_adjs2, UNION);
        std::sort(tmp_adjs2.begin(), tmp_adjs2.end());
        tmp_adjs.resize(tmp_adjs2.size());
        tmp_adjs.erase(std::set_intersection(b_ents.begin(), b_ents.end(),
                                             tmp_adjs2.begin(), tmp_adjs2.end(),
                                             tmp_adjs.begin()), tmp_adjs.end());
      }

      if (2 == dimension() && tmp_adjs.size() > 1) {
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

void ModelEnt::get_indexed_connect_coords(std::vector<moab::EntityHandle> &ents,
                                          std::vector<int> *senses,
                                          moab::Tag tagh,
                                          std::vector<int> &ents_as_ids,
                                          std::vector<double> &coords,
                                          moab::Range *verts_range) 
{
    // if we need to, make a tag
  moab::ErrorCode rval;
  bool i_created_tag = false;
  if (0 == tagh) {
    rval = mk_core()->moab_instance()->tag_create("__ModelEntidtag", sizeof(int), moab::MB_TAG_DENSE, 
                                                  moab::MB_TYPE_INTEGER, tagh, NULL);
    MBERRCHK(rval, mkCore->moab_instance());
    i_created_tag = true;
  }
  
    // put entities into range, after clearing it
  moab::Range tmp_range, ents_range;
  if (!verts_range) verts_range = &tmp_range;
  verts_range->clear();

  std::copy(ents.begin(), ents.end(), moab::range_inserter(ents_range));
  bool all_verts = (mkCore->moab_instance()->type_from_handle(*ents_range.begin()) ==
                    mkCore->moab_instance()->type_from_handle(*ents_range.begin()) &&
                    mkCore->moab_instance()->type_from_handle(*ents_range.begin()) == moab::MBVERTEX);


    // get connectivity of all ents and store in range
  if (!all_verts) {
    rval = mkCore->moab_instance()->get_connectivity(ents_range, *verts_range);
    MBERRCHK(rval, mkCore->moab_instance());
  }
  else *verts_range = ents_range;

    // resize index array, to max number of connectivity entries
  int max_numconnect = moab::CN::VerticesPerEntity(mkCore->moab_instance()->type_from_handle(*(ents_range.rbegin())));
  ents_as_ids.resize(ents.size()*max_numconnect);
  
    // temporarily store ids in ents_as_ids array, and set id tag 
  assert(ents_as_ids.size() >= verts_range->size());
  for (unsigned int i = 0; i < verts_range->size(); i++) ents_as_ids[i] = i;
  rval = mk_core()->moab_instance()->tag_set_data(tagh, *verts_range, &ents_as_ids[0]);
  MBERRCHK(rval, mkCore->moab_instance());

    // get ids into ids vector in connectivity or ents order
  if (all_verts) {
    rval = mk_core()->moab_instance()->tag_get_data(tagh, &ents[0], ents.size(), &ents_as_ids[0]);
    MBERRCHK(rval, mkCore->moab_instance());
    ents_as_ids.resize(ents.size());
  }
  else {
      // loop over entities
    std::vector<moab::EntityHandle> storage;
    const moab::EntityHandle *connect;
    int num_connect;
    unsigned int last = 0;
    for (unsigned int i = 0; i < ents.size(); i++) {
        // first, get connect vector
      mk_core()->moab_instance()->get_connectivity(ents[i], connect, num_connect, true, &storage);
      MBERRCHK(rval, mk_core()->moab_instance());
      assert(ents_as_ids.size() >= last + num_connect);
        // next, get ids and put into ids vector
      rval = mk_core()->moab_instance()->tag_get_data(tagh, connect, num_connect, &ents_as_ids[last]);
      MBERRCHK(rval, mk_core()->moab_instance());
        // reverse if necessary
      assert(!senses || (*senses)[i] != SENSE_BOTH);
      if (senses && (*senses)[i] == SENSE_REVERSE) std::reverse(&ents_as_ids[last], &ents_as_ids[last+num_connect]);
        // update last
      last += num_connect;
    }
  }

    // get coords for range-ordered vertices
  coords.resize(3*verts_range->size());
  rval = mk_core()->moab_instance()->get_coords(*verts_range, &coords[0]);
  MBERRCHK(rval, mk_core()->moab_instance());

    // delete the tag if I created it
  if (i_created_tag) {
    rval = mk_core()->moab_instance()->tag_delete(tagh);
    MBERRCHK(rval, mkCore->moab_instance());
  }
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
  iRel::Error err = mkCore->irel_pair()->getEntSetRelation(gent, false, IBSHR(ment));
  IBERRCHK(err, "Failed to get mesh set handle for geometry entity.");
  return ment;
}

    /** \brief Get mesh interval size, if any
     * Returns -1 if no size set on this entity.  If intervals are set, returns computed size.
     * \return Interval size for this ModelEnt.
     */
double ModelEnt::mesh_interval_size() const 
{
  if (-1 != sizing_function_index() && mk_core()->sizing_function(sizing_function_index()) &&
      -1 != mk_core()->sizing_function(sizing_function_index())->size())
    return mk_core()->sizing_function(sizing_function_index())->size();
  
  else if (1 == dimension() && -1 != mesh_intervals() && 0 != mesh_intervals())
    return measure() / mesh_intervals();
    
  else return -1;
}

} // namespace MeshKit

