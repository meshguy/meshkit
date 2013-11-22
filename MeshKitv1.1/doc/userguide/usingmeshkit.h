/*!
  \page usingmeshkit Using MeshKit: The 2-Minute Introduction
  In its simplest form, %MeshKit can be used by specifying the geometry, a mesh scheme,
  and a desired mesh size, then executing the meshing operation.  In code, this looks like:
  \code
  MeshKit::MKCore mk;                             // by default, creates geometry, mesh, and relations instances automatically
  mk.load_geometry(filename);                     // load a geometric model into MeshKit
  MeshKit::MEntVector model_ents;
  mk.get_entities_by_dimension(3, vols);          // get all geometric volumes
  MeshKit::SizingFunction esize(mk, -1, 0.25);    // create a sizing function to use everywhere
  mk.construct_meshop("NGTetMesher", vols);       // construct a tet mesher and put it in the graph
  mk.setup_and_execute();                         // execute the meshing graph
  \endcode
  The mesh operation, here the Netgen tet mesher, understands how to "setup" or prepare the volume
  for tet meshing.  Behind the scenes, a tri mesher is added to the graph, then setup is called on that
  tri mesher; it in turn creates an edge mesher and assigns all geometric edges to it.  If not specified
  directly, the size used to mesh each entity is the first one registered with the MeshKit instance.
  Once generated, the mesh can be retrieved with various functions on the MeshKit instance, or 
  directly through the MOAB or iMesh interfaces.

  For more complicated meshing processes, a directed graph of mesh-based operations can be specified, with
  each graph node representing an operation on one or more model entities.  This can represent
  a simple topology-driven meshing process, where entities are meshed in order of increasing 
  topological dimension; or it can be a sequence of more general steps, for example a mesh, smooth,
  refine sequence.  The graph can be set up before any meshing is done.  Nodes can be invalidated 
  (e.g. by deleting mesh produced by that node), which invalidates later nodes in the graph, and 
  the graph can be re-executed.  For more information on the graph-based meshing process implemented
  in MeshKit, see \ref thegraph.
 
Top: \ref index Prev: \ref interfaces Next: \ref api
 */
