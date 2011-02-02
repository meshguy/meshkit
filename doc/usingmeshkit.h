/*!
  \page usingmeshkit Using %MeshKit
  In its simplest form, %MeshKit can be used simply by specifying the geometry, a mesh scheme,
  and a desired mesh size, then executing the meshing operation.  In code, this looks like:
  \code
  MeshKit::MKCore mk; // by default, creates geometry, mesh, and relations instances automatically
  mk.load_geometry(filename);
  MeshKit::MEVector model_ents;
  mk.get_entities_by_dimension(3, model_ents);
  MeshKit::TetMeshScheme tm(model_ents, SizingFunction(0.25));
  tm.execute();
  \endcode
  For more complicated meshing processes, a directed graph of meshing operations can be specified, with
  each graph node representing a meshing operation on one or more model entities.  The graph can
  be set up before any meshing is done.  Nodes can be invalidated (e.g. by deleting mesh produced
  by that node), which invalidates later nodes in the graph, and the graph can be re-executed.
 
Top: \ref index Prev: \ref interfaces Next: \ref api
 */
