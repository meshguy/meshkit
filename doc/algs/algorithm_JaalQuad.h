/*!
 \page algorithm_JaalQuad Jaal Quad Mesher

<b>Name:</b> QuadMesher

<b>External dependencies:</b> Mesquite (optional)

<b>Input:</b> 2D ModelEnt's, <em>meshed with triangles</em>

<b>Output:</b> Mesh vertices, quadrilaterals, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

Looks on the ModelEnt for an assigned tri mesher, and if none is found, assigns one from the default list
of tri meshing algorithms registered with MeshKit.  <em>No check is made on the number of bounding intervals
being even; if this number is odd, one triangle will remain after this algorithm is finished.</em>

<b>Notes:</b>

This class implements quadrilateral mesh generation, based on triangle recombination.  Triangles can be
generated with any algorithm.  Irregular nodes are removed using Bunin's algorithm.

*/
