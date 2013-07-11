/*!
  \page interfaces Geometry, Mesh, Relations Interfaces
  %MeshKit interacts with geometry, mesh, and relations through well-defined APIs for each.  A %MeshKit
  instance keeps a reference to a single instance of each of those APIs, and stores all generated
  data in those interface instances.  For geometry and relations, %MeshKit uses the ITAPS %iGeom and %iRel
  interfaces; for mesh, portions of %MeshKit use both the MOAB and ITAPS %iMesh interfaces.  The C-based
  ITAPS interfaces can be rather cumbersome to use, especially when accessing lists of entities (a quite 
  common operation in mesh generation).  Therefore, %MeshKit defines the iGeom, iMesh, and iRel C++ classes 
  for the ITAPS interfaces; these classes are strictly wrappers, and call the C-based functions directly.  
  Lists are converted to/from std::vectors, eliminating the need for careful memory management of these 
  lists in the mesh classes.
 
  At first glance, the choice to use the MOAB interface in %MeshKit would seem to prevent interoperability
  with other mesh databases, contradicting one of the primary goals of the %MeshKit design.
  The are two primary reasons for this choice.  First, the use of MOAB's Range class
  can result in substantial savings, in both memory and execution time (see the MOAB User's Guide for
  more discussion of this class).  Second, MOAB provides important utilities not available
  through %iMesh, and which are used in several important parts of %MeshKit; examples include a skinning
  tool, various types of tree decompositions, and utilities for accessing geometric topology
  stored in the mesh entity sets.
 
  To preserve interoperability, then, we intend to implement MOAB's interface on top of %iMesh.  While this
  may result in sub-optimal performance on other databases, we hope the overhead will not prove too severe.

Top: \ref index Prev: \ref thegraph Next: \ref usingmeshkit
  
 */
