/** \mainpage %MeshKit: A Library for Mesh Generation
 *
 * %MeshKit is a library for geometry-based 2d and 3d mesh generation.  This library has two main purposes:
 * - to provide a useful collection of open-source mesh generation functions
 * - to provide infrastructure (geometry, lower-dimensional entity meshing, etc.) for implementing new
 *   meshing functionality
 *
 * %MeshKit uses a component-based approach, using external components for selected functionality
 * where possible.  Geometry and relations to mesh are accessed through the %iGeom and %iRel interfaces,
 * resp.  Mesh is accessed through both the MOAB and %iMesh interfaces.  %MeshKit allows registration
 * of external meshing algorithms, allowing those algorithms to operate in collaboration with the
 * rest of %MeshKit functionality.
 *
 * \section datamodel Data Model
 * %MeshKit defines the following top-level classes:
 * - <b> MKCore: </b> Top-level instance providing access to most other data in %MeshKit
 * - <b> ModelEntity: </b> Geometric model entity and the mesh associated with it
 * - <b> MeshOp: </b> An operation associated with the mesh generation process
 * - <b> MeshScheme: </b> A type of MeshOp that generates mesh
 * - <b> SizingFunction: </b> An object storing information about mesh size specification
 *
 * The data model used by %MeshKit is similar to other geometry-based mesh generation tools, revolving 
 * around geometric model entities and the algorithms used to mesh them.  Data is stored down in
 * the geometry and mesh interface implementations, then referenced from objects in %MeshKit.
 * 
 * \section interfaces Geometry, Mesh, Relations Interfaces
 * %MeshKit interacts with geometry, mesh, and relations through well-defined APIs for each.  A %MeshKit
 * instance keeps a reference to a single instance of each of those APIs, and stores all generated
 * data in those interface instances.  For geometry and relations, %MeshKit uses the ITAPS %iGeom and %iRel
 * interfaces; for mesh, portions of %MeshKit use both the MOAB and ITAPS %iMesh interfaces.  The C-based
 * ITAPS interfaces can be rather cumbersome to use, especially when accessing lists of entities (a quite 
 * common operation in mesh generation).  Therefore, %MeshKit defines the iGeom, iMesh, and iRel C++ classes 
 * for the ITAPS interfaces; these classes are strictly wrappers, and call the C-based functions directly.  
 * Lists are converted to/from std::vectors, eliminating the need for careful memory management of these 
 * lists in the mesh classes.
 *
 * At first glance, the choice to use the MOAB interface in %MeshKit would seem to prevent interoperability
 * with other mesh databases, contradicting one of the primary goals of the %MeshKit design.
 * The are two primary reasons for this choice.  First, the use of MOAB's Range class
 * can result in substantial savings, in both memory and execution time (see the MOAB User's Guide for
 * more discussion of this class).  Second, MOAB provides important utilities not available
 * through %iMesh, and which are used in several important parts of %MeshKit; examples include a skinning
 * tool, various types of tree decompositions, and utilities for accessing geometric topology
 * stored in the mesh entity sets.
 *
 * To preserve interoperability, then, we intend to implement MOAB's interface on top of %iMesh.  While this
 * may result in sub-optimal performance on other databases, we hope the overhead will not prove too severe.
 * 
 * \section usingmeshkit Using MeshKit
 * In its simplest form, %MeshKit can be used simply by specifying the geometry, a mesh scheme,
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
 * %MeshKit is implemented in C++, but also includes auto-generated python interfaces to all API 
 * classes.  Access is through the MKCore class, which is instantiated by the application.  All %MeshKit
 * classes are declared inside the MeshKit namespace.
 *
 */
