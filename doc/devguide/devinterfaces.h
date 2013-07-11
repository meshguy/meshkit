/*!
  \page devinterfaces Geometry, Mesh, Relations Interfaces
The MKCore instance stores pointers to the geometry, mesh, and relations interfaces used elsewhere in MeshKit,
as well as member functions for accessing those instance pointers.  Individual algorithms do not need to store
their own pointers or references to those instances.  Classes derived from MeshOp also do not need to store a 
pointer to MKCore, since one is already stored in GraphNode.

The following member functions should be used to get various interface instances.
 - MKCore* MeshOp::mk_core(): Returns a MKCore*; this is a true function call, and performs a dyanamic_cast
    from a MKGraph object, which is actually stored in GraphNode.
 - iGeom* MKCore::igeom_instance(int index=0)
 - iGeom_Instance iGeom::instance()
 - iMesh* MKCore::imesh_instance(int index=0)
 - iMesh_Instance iMesh::instance()
 - moab::Interface* MKCore::moab_instance(int index=0)
 - iRel* MKCore::irel_instance(int index=0)
 - iRel_Instance iRel::instance()

Note that MKCore can store pointers to multiple iGeom, iMesh, or iRel instances.  Most of the time, though, a single 
instance of these interfaces is used.

Top: \ref index Prev: \ref api Next: \ref styleguide
  
 */
