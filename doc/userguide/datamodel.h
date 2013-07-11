/*!
  \page datamodel Data Model
  %MeshKit defines the following top-level classes:
  - <b> MKCore: </b> Top-level instance providing access to most other data in %MeshKit
  - <b> ModelEntity: </b> Geometric model entity and the mesh associated with it
  - <b> MeshOp: </b> An operation associated with the mesh generation process
  - <b> MeshScheme: </b> A type of MeshOp that generates mesh
  - <b> SizingFunction: </b> An object storing information about mesh size specification
 
  The data model used by %MeshKit is similar to other geometry-based mesh generation tools, revolving 
  around geometric model entities and the algorithms used to mesh them.  Data is stored down in
  the geometry and mesh interface implementations, then referenced from objects in %MeshKit.

Top: \ref index Prev: \ref index Next: \ref thegraph "A Graph-Based Model for Mesh Generation"
  
 */
