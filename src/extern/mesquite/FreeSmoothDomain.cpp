#include "meshkit/ModelEnt.hpp"
#include "meshkit/MKCore.hpp"
#include "moab/Interface.hpp"
#include "meshkit/FreeSmoothDomain.hpp"
#include "moab/mesquite/MsqError.hpp"

namespace MeshKit {

using namespace moab;
using namespace MESQUITE_NS;

FreeSmoothDomain::FreeSmoothDomain( MKCore* core, const MEntVector& entities )
  : MsqCommonIGeom(core->igeom_instance()->instance()),
    haveEntGeomRelTag(false),
    moabIface(core->moab_instance())
{
    // Group input entities and all child entities by dimension.
    // This way if a vertex is contained in all entities in which it
    // is used rather than just the one that owns it, it will still
    // get assigned correctly as long as we work from highest to lowest
    // dimension.
  MEntVector ents_by_dim[4], tmp_vec;
  MEntVector::const_iterator i;
  for (i = entities.begin(); i != entities.end(); ++i) {
    ModelEnt* ent = *i;
    int dim = ent->dimension();
    if (dim < 0 || dim > 3)
      throw Error(MK_WRONG_DIMENSION, "Entity of invalid dimension %d", dim);
    ents_by_dim[dim].push_back(ent);
    for (int d = 0; d < dim; ++d) {
      tmp_vec.clear();
      ent->get_adjacencies( d, tmp_vec );
      ents_by_dim[d].insert( ents_by_dim[d].end(), tmp_vec.begin(), tmp_vec.end() );
    }
  }
  
    // Create a tag to store geometry handles on mesh entities
  assert( sizeof(iBase_EntityHandle) == sizeof(moab::EntityHandle) );
  const moab::EntityHandle zero = 0;
  moab::ErrorCode rval = moabIface->tag_get_handle( 0, 
                                          1,
                                          MB_TYPE_HANDLE,
                                          entGeomRel,
                                          moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
                                          &zero );
  MBERRCHK( rval, moabIface );
  haveEntGeomRelTag = true;

    // NOTE: We only go up to a dimension of 2 here!
    //       Mesquite only cares about surface elements constrained
    //       to a surface and vertices constrained to lie on surfaces,
    //       curves, or points.
  for (int d = 0; d < 3; ++d) {
    std::sort( ents_by_dim[d].begin(), ents_by_dim[d].end() );
    MEntVector::iterator end = std::unique( ents_by_dim[d].begin(), ents_by_dim[d].end() );
    for (i = ents_by_dim[d].begin(); i != end; ++i) {
      EntityHandle set = (*i)->mesh_handle();
      iBase_EntityHandle geom = (*i)->geom_handle();
      if (!geom) 
        throw Error(MK_BAD_GEOMETRIC_EVALUATION, "Cannot do free smooth if ModelEnt doesn't have geometry" );
      
        // assign surface elements
      if (d == 2) { 
        Range elems;
        rval = moabIface->get_entities_by_dimension( set, 2, elems );
        MBERRCHK( rval, moabIface );
        rval = moabIface->tag_clear_data( entGeomRel, elems, &geom );
        MBERRCHK( rval, moabIface );
      }
      
        // assign vertices
      Range verts;
      rval = moabIface->get_entities_by_dimension( set, 0, verts );
      MBERRCHK( rval, moabIface );
      rval = moabIface->tag_clear_data( entGeomRel, verts, &geom );
      MBERRCHK( rval, moabIface );
    }
  }
}

FreeSmoothDomain::~FreeSmoothDomain()
{
  if (haveEntGeomRelTag) 
    moabIface->tag_delete( entGeomRel );
}

iBase_EntityHandle FreeSmoothDomain::get_geometry( Mesh::EntityHandle mesh_ent ) const
{
  iBase_EntityHandle result;
  moab::EntityHandle ent = (moab::EntityHandle)mesh_ent;
  moab::ErrorCode rval = moabIface->tag_get_data( entGeomRel, &ent, 1, &result );
  MBERRCHK( rval, moabIface );
  return result;
}

void FreeSmoothDomain::get_geometry( const Mesh::EntityHandle* mesh_ents,
                                     size_t num_handles,
                                     iBase_EntityHandle* geom_ents,
                                     MsqError& err ) const
{
  const moab::EntityHandle* ents = (const moab::EntityHandle*)mesh_ents;
  moab::ErrorCode rval = moabIface->tag_get_data( entGeomRel, ents, num_handles, geom_ents );
  if (moab::MB_SUCCESS != rval) {
    MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,
                    "Error querying MOAB for geom");
  }                    
}


void FreeSmoothDomain::snap_to( Mesh::VertexHandle handle,
                                Vector3D& coordinate ) const
{
  iBase_EntityHandle geom = get_geometry( handle );
  if (geom) {
    int ierr = move_to( geom, coordinate );
    IBERRCHK(ierr, "iGeom evaluation error");
  }
}

void FreeSmoothDomain::vertex_normal_at( Mesh::VertexHandle handle,
                                         Vector3D& coordinate ) const
{
  iBase_EntityHandle geom = get_geometry( handle );
  if (geom) {
    int ierr = normal( geom, coordinate );
    IBERRCHK(ierr, "iGeom evaluation error");
  }
}

void FreeSmoothDomain::element_normal_at( Mesh::ElementHandle handle,
                                     Vector3D& coordinate ) const
{
  FreeSmoothDomain::vertex_normal_at( handle, coordinate );
}

void FreeSmoothDomain::vertex_normal_at( const Mesh::VertexHandle* handle,
                                         Vector3D coordinates[],
                                         unsigned count,
                                         MsqError& err ) const
{
  tmpHandles.resize( count );
  get_geometry( handle, count, &tmpHandles[0], err );
  MSQ_ERRRTN(err);
  
  int ierr = normal( &tmpHandles[0], coordinates, count );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,"Failure to evaluate geometry.  iBase ErrorCode %d\n", ierr);
}

void FreeSmoothDomain::domain_DoF( const Mesh::VertexHandle* handle_array,
                                   unsigned short* dof_array,
                                   size_t count,
                                   MsqError& err ) const
{
  tmpHandles.resize( count );
  get_geometry( handle_array, count, &tmpHandles[0], err );
  MSQ_ERRRTN(err);

  int ierr, type;
  for (size_t i = 0; i < count; ++i) {
    if (!tmpHandles[i])
      dof_array[i] = 3; // don't have geom handles for volumes
    else if (i > 0 && tmpHandles[i-1] == tmpHandles[i])
      dof_array[i] = dof_array[i-1];
    else {
      iGeom_getEntType( geomIFace, tmpHandles[i], &type, &ierr );
      dof_array[i] = type;
      
      if (iBase_SUCCESS != ierr) {
        MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,"Failure to evaluate geometry.  iBase ErrorCode %d\n", ierr);
        return;
      }
    }
  }
}

void FreeSmoothDomain::closest_point( Mesh::VertexHandle handle,
                                      const Vector3D& position,
                                      Vector3D& closest,
                                      Vector3D& normal,
                                      MsqError& err ) const
{
  iBase_EntityHandle geom = get_geometry( handle );
  if (geom) {
    int ierr = closest_and_normal( geom, position, closest, normal );
    if (iBase_SUCCESS != ierr) {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR,"Failure to evaluate geometry.  iBase ErrorCode %d\n", ierr);
    }
  }
}


} // namespace MeshKit
