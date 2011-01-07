/** \mainpage MeshKit: A Library for Mesh Generation
 *
 * MeshKit is a library for geometry-based 2d and 3d mesh generation.  This library has two main purposes:
 * - to provide a useful collection of open-source mesh generation functions
 * - to provide infrastructure (geometry, lower-dimensional entity meshing, etc.) for implementing new
 *   meshing functionality
 *
 * MeshKit uses a component-based approach, using external components for selected functionality
 * where possible.  Geometry and relations to mesh are accessed through the iGeom and iRel interfaces,
 * resp.  Mesh is accessed through the MOAB interface; an implementation of that interface on top of
 * iMesh is provided, for those wanting to use a different mesh database.  MeshKit allows registration
 * of external meshing algorithms, allowing those algorithms to operate in collaboration with the
 * rest of MeshKit functionality.
 *
 * \section datamodel Data Model
 * MeshKit defines the following top-level classes:
 * - <b> MKCore: </b> Top-level instance providing access to most other data in MeshKit
 * - <b> ModelEntity: </b> Geometric model entity and the mesh associated with it
 * - <b> MeshOp: </b> An operation associated with the mesh generation process
 * - <b> MeshScheme: </b> A type of MeshOp that generates mesh
 * - <b> SizingFunction: </b> An object storing information about mesh size specification
 *
 * The data model used by MeshKit is similar to other geometry-based mesh generation tools, revolving 
 * around geometric model entities and the algorithms used to mesh them.  Data is stored down in
 * the geometry and mesh interface implementations, then referenced from objects in MeshKit.
 * 
 * \section usingmeshkit Using MeshKit
 * In its simplest form, MeshKit can be used simply by specifying the geometry, a mesh scheme,
 * and a desired mesh size, then executing the meshing operation.  In code, this looks like:
 * \code
 * MeshKit::MKCore mk; // by default, creates geometry, mesh, and relations instances automatically
 * mk.load_geometry(filename);
 * MeshKit::MEVector model_ents;
 * mk.get_entities_by_dimension(3, model_ents);
 * MeshKit::TetMeshScheme tm(model_ents, SizingFunction(0.25));
 * tm.execute();
 * \endcode
 * For more complicated meshing processes, a directed graph of meshing operations can be specified, with
 * each graph node representing a meshing operation on one or more model entities.  The graph can
 * be set up before any meshing is done.  Nodes can be invalidated (e.g. by deleting mesh produced
 * by that node), which invalidates later nodes in the graph, and the graph can be re-executed.
 *
 * \section api MeshKit API
 * MeshKit is implemented in C++, but also includes auto-generated python interfaces to all API 
 * classes.
 *
 */
