/*! \mainpage MeshKit: A Library for Mesh Generation

MeshKit is a library for geometry-based 2d and 3d mesh generation.  This library has two main purposes:
 - to provide a useful collection of open-source mesh generation functions
 - to provide infrastructure (geometry, lower-dimensional entity meshing, etc.) for implementing new
   meshing functionality

 %MeshKit uses a component-based approach, using external components for selected functionality
 where possible.  Geometry and relations to mesh are accessed through the %iGeom and %iRel interfaces,
 resp.  Mesh is accessed through both the MOAB and %iMesh interfaces.  %MeshKit allows registration
 of external meshing algorithms, allowing those algorithms to operate in collaboration with the
 rest of %MeshKit functionality.

 This manual is divided into the following sections:

- \subpage userguide
  - \ref usingmeshkit
  - \ref datamodel
  - \ref thegraph
    - \ref detailedtraversal
    - \ref graphexamples
  - \ref interfaces
  - \ref api
  - \ref meshkitalgorithms
  - <a href="examples.html"><b>Examples</a>

- \subpage devguide
  - \ref styleguide
  - \ref devinterfaces
  - \ref newmeshop
  - \ref faq

 */
