/* \page copyexpandextrude_sets Copy, Expand, and Extrude Sets

When copying or extruding mesh, it is often desired to also handle boundary condition or other
groupings associated with the mesh.  MeshKit handles this through definition of Copy, Expand, and Extrude sets.
Mesh contained in each of these types of sets is handled as follows:

 - <b>Copy sets:</b> For each copy set A, a new copy set A' is created.  For each entity e (being copied to e')
   contained in A, the copy e' is put in A'.  An optional tag points from A to A'.

 - <b>Expand sets:</b> For each expand set A, if an entity e (being copied to e') is in A, its copy e' is also 
   put into A.

 - <b>Extrude sets:</b> For each extrude set A, containing an entity e (being extruded to a set of entities {e'},
   dim(e') = dim(e)+1), e is \e replaced in A by e'.  Thus, if set A contains entities with dimension d before
   the extrude, it ends up with entities with dimension d+1.  The number of members of A after the extrude is
   N-1 times the number of members before the extrude, where N is the number of extrude steps.

Copy, expand, and extrude sets enable mesh generation algorithms to be implemented in a way that is independent
of the semantics of a given type of grouping.  For example, material specifications, boundary condition grouping,
and geometric model relationships can all be updated through copy, extrude, or refinement operations, without
those operations understanding the notion of materials, boundary conditions, or geometric modeling.

*/
