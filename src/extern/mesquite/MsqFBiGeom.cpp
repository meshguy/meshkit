#include "meshkit/MKCore.hpp"
#include "meshkit/MsqFBiGeom.hpp"

namespace MeshKit
{

MsqFBiGeom::MsqFBiGeom( FBiGeom * ifbigeom,
            iBase_EntityHandle geom_ent_handle )
{
  geomEntHandle = geom_ent_handle;
  fbigeom = ifbigeom;
}

MsqFBiGeom::~MsqFBiGeom() {}

void MsqFBiGeom::snap_to( Mesh::VertexHandle ,
    Vector3D& coord ) const
{
  double x, y, z;
  iBase_ErrorType err = fbigeom->getEntClosestPt(geomEntHandle, coord[0], coord[1], coord[2], x, y, z );
  IBERRCHK(err, "Failed to get closest point to surface.");
  coord.set( x, y, z );
}

void MsqFBiGeom::vertex_normal_at( Mesh::VertexHandle /*entity_handle*/,  // the mesh vertex handle is unused
                         Vector3D& coord ) const
{
  double i, j, k;
  iBase_ErrorType err = fbigeom->getEntNrmlXYZ( geomEntHandle, coord[0], coord[1], coord[2], i, j, k);
  IBERRCHK(err, "Failed to get normal to surface.");
  coord.set( i, j, k );

}

void MsqFBiGeom::element_normal_at( Mesh::ElementHandle /*entity_handle*/,  // the mesh element handle is unused
                          Vector3D& coord ) const
{
  double i, j, k;
  iBase_ErrorType err = fbigeom->getEntNrmlXYZ( geomEntHandle, coord[0], coord[1], coord[2], i, j, k);
  IBERRCHK(err, "Failed to get normal to surface.");
  coord.set( i, j, k );
}

void MsqFBiGeom::vertex_normal_at( const Mesh::VertexHandle* /*handles*/,
                         Vector3D coordinates[],
                         unsigned count,
                         MsqError& err ) const
{
  for (unsigned int i=0; i<count; i++)
  {
    double x, y, z;
    Vector3D coord =coordinates[i];
    iBase_ErrorType err = fbigeom->getEntNrmlXYZ( geomEntHandle, coord[0], coord[1], coord[2], x, y, z);
    IBERRCHK(err, "Failed to get normal to surface.");
    coordinates[i].set( x, y, z );
  }

}

void MsqFBiGeom::closest_point( Mesh::VertexHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const
{
  iBase_ErrorType ierr = fbigeom->getEntNrmlPlXYZ( geomEntHandle,
                           position[0], position[1], position[2],
                           closest[0], closest[1], closest[2],
                           normal[0],  normal[1],  normal[2]);
  IBERRCHK(ierr, "Failed to get normal and closest point to surface.");
}

void MsqFBiGeom::domain_DoF( const Mesh::VertexHandle* handle_array,
                   unsigned short* dof_array,
                   size_t num_vertices,
                   MsqError& err ) const
{
  // just fill with 2, this is a surface, we know
  for (size_t i = 0; i< num_vertices; i++)
  {
    dof_array[i]=2;
  }
}

}// end namespace
