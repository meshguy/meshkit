namespace MeshKit {

/*!
  \page newmeshop Writing A New MeshOp Class
%MeshKit is designed to simplify development of new meshing algorithms, by providing support for
 - Gathering boundary mesh from lower-dimensional bounding geometry,
 - Coordinating relations between geometry and mesh,
 - Setup and execution of a graph-based meshing process,

and many other things one would otherwise need to do to support their new algorithm.  In this section, we describe the
basic steps required to write a new meshing algorithm and integrate it with the rest of %MeshKit.

\section overview Overview
In %MeshKit, each meshing operation is implemented in a class derived from MeshOp, which provides support for executing
the new algorithm in a graph-based meshing process.  The meshing operation class is registered with MKCore, which allows 
the algorithm to be requested by name in applications.  The operation can provide a function which evaluates whether
a given entity can be meshed with that algorithm, providing further automation of the meshing process.  Each meshing operation 
instance stores a map between ModelEnt objects and moab::Range's operated on by the operation; this can be used as input 
to the next node in the meshing graph, or to reverse the effect of the MeshOp.

The following sections describe functions that must be implemented for a new MeshOp in order for it to integrate into %MeshKit.

\section registration Declaration & Registration
Each new meshing operation should be derived from the MeshOp class (see e.g. VertexMesher for the proper syntax).  The new MeshOp
should be registered with MKCore at some point before an application requests the new MeshOp.  For MeshOp classes packaged
with MeshKit, this happens inside a few static functions called during MeshKit startup.  For application-provided MeshOps,
registration can happen either during static initializaton of the application, or early in the application's execution
before the new MeshOp is needed.  

The following steps should be taken to register a new MeshOp with MKCore:
 -# For each pure virtual function defined in MeshOpProxy, a corresponding
   static function should be defined in the new MeshOp class, with the
   exception of the \c create method, for which the MeshOp should provide
   a constructor that accepts the same argument list.  The list of such
   static functions is currently:
   - <tt> const char* name() </tt> - The meshing operation name (used by applications to request this meshing operation).
   - <tt> bool can_mesh(iBase_EntityType dimension) </tt> - Returns true
      for any topogical dimension(s) of ModelEnt's this meshing operation 
      can operate on. <em> The behavior of this overloaded variant of <tt> can_mesh </tt> is different from the version
      taking a ModelEnt* as an argument; see \ref canmesh below.</em>
   - <tt> const moab::EntityType* output_types() </tt> - Returns an array
      of mesh entity types terminated with a value of  moab::MBMAXTYPE .  
      This function does not return the size of the array, so the terminating
      value must be present for calling code to know where the end of the
      array is.  The values contained in this array (with the exception of the
      terminating value) are the types of mesh entities produced or operated on by the meshing
      algorithm.
   - <tt> bool can_mesh(ModelEnt* ent) </tt> - evaluates whether the operation can operate on a specified ModelEnt.
      The MeshOp class provides several convenience functions that can be used for the canmesh function.  These functions
      are canmesh_vertex, canmesh_edge, canmesh_surface, canmesh_region, and return true if the specified ModelEnt has topological
      dimension 0, 1, 2, or 3, respectively.
 -# Add a line to either algs/register_algs.cpp or extern/register_extern.cpp
   of the form <tt> REGISTER_MESH_OP(classname); </tt>.  If adding a
   MeshOp in \c extern/, use preprocesor directives as appropriate to ensure
   that %MeshKit builds correctly if configured without support for the 
   corresponding external libraries.

Classes derived from MeshOp should also provide the mesh_types_arr function used by the MeshOp base class to 
implement the mesh_types funciton.  This can typically be implemented as:
\code
  virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }
\endcode

If for some reason the implementation of static methods as decribed above
or the RegisterMeshOp template used by the REGISTER_MESH_OP(classname) macro are inappropriate
for a given MeshOp, you may alternately implement your own subclass of 
MeshOpProxy for your MeshOp and explcitly register by calling MeshKit::register_meshop
from either \c register_algs_mesh_ops() defined in \c algs/register_algs.cpp
or \c register_extern_mesh_ops() defined in \c extern/register_extern.cpp.  Alternatively, an application can register
new MeshOps early in application execution, using the same macros used in \c register_algs_mesh_ops().

\section canmesh can_mesh Functions

Two variants of the \c can_mesh function are used in MeshOp registration.  The first variant, 
 <tt> bool can_mesh(iBase_EntityType dimension) </tt>, returns whether it is <em> possible</em> to mesh ModelEnt's with the specified
dimension with the MeshOp.  This evaulation is made purely on the basis of topological dimension.  In some cases, there
may be further constraints on the MeshOp such that, given a concrete ModelEnt, the meshOp will not be able to mesh the
ModelEnt.  For example, a surface mapping algorithm would return \c true for iBase_EntityType=iBase_FACE, but when passed
a ModelEnt corresponding to a non-simply-connected surface (a surface with holes) would return \c false.

The MeshOp class provides several convenience functions that can be used in place of a MeshOp-specific \c can_mesh function.  
These functions are \c canmesh_vertex, \c canmesh_edge, \c canmesh_surface, \c canmesh_region, and return \c true if the 
specified \c ModelEnt has topological dimension 0, 1, 2, or 3, respectively.  These functions can be passed directly to the 
registration code, so the new MeshOp need not implement its own \c can_mesh function.

\section setupexecute setup_this and execute_this
Two functions, MeshOp::setup_this and MeshOp::execute_this, must be implemented to support the graph-based meshing process in %MeshKit.

setup_this performs setup-type tasks for the associated ModelEnt(s), like setting intervals or computing vertex or edge types.  
When applied to geometric model entities, it is also responsible for coordinating setup of lower-dimensional bounding entities.
It does this by placing more meshing operations on the graph, upstream of the current MeshOp.  For example, a surface meshing
MeshOp would put MeshOp nodes on the graph for meshing the surface's bounding edges; the MeshOps for the bounding edges
would create MeshOp nodes for meshing the bounding vertices.  As a general rule, though, setup_this
should not create any mesh entities.

In many cases, %MeshKit can automatically choose an appropriate MeshOp for meshing bounding entities, based on
the canmesh functions specified at registration.  The MeshOp::setup_boundary() function performs this task, placing resulting
graph nodes on the graph upstream of the current MeshOp.  setup_boundary() can be called
directly by other MeshOps.  For example, the EdgeMesher class in %MeshKit calls this function to accomplish the trivial task
of setting up meshing for geometric vertices bounding the edge.

The execute_this() function should perform the actual meshing operation.  This function is called starting at the root of the
graph and proceeding in the forward direction.  Most often, execute_this() is called on an entity only after it has been 
called for all its lower-dimensional bounding entities.

After the mesh is generated for an entity, it should be "committed" to the mesh database and associated with the geometric
model entity.  This is most easily accomplished by calling the MeshOp::commit_mesh function, passing only the entities
created or operated on by the MeshOp.  For example, in EdgeMesh, only the mesh vertices and edges owned by the model edge
are passed to commit_mesh; this mesh does not include the vertices bounding the model edge.  execute_this marks the ModelEnt 
with a MeshedState corresponding to unmeshed, partially meshed, or fully meshed state.

\section othersupport Other Support or Convenience Functions

Many other functions are provided to simplify development of meshing operations:
 - ModelEnt::boundary() gathers bounding mesh or model entities, returning them in various types of lists;
 this function comes in various overloaded forms, depending on the type of list returned.
 - ModelEnt::get_adjacencies provides adjacency lists in various forms.
 - Geometric evaluations are provided by ModelEnt::measure and ModelEnt::evaluate.  The scope of functions like this is
 intentionally narrow, to avoid duplication of the full geometric model interface available in iGeom.  iGeom can be used
 to perform other types of geometric or topoligical evaluations not directly available from ModelEnt.
 

Top: \ref index Prev: \ref api Next: \ref styleguide
  
 */

}
