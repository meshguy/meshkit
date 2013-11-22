/*!
 \page algorithm_edgemesher EdgeMesher

<b>Name:</b> EdgeMesher

<b>External dependencies:</b> (none)

<b>Input:</b> 1D ModelEnt's

<b>Output:</b> Mesh vertices, edges, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

Calls MeshOp::setup_boundary().

<b>Notes:</b>

This class implements four specific edge meshing schemes:

 - EQUAL: equal spacing along arc length
 - BIAS: bias spacing toward either end
 - DUALBIAS: bias spacing towards both ends
 - CURVATURE: curvature-based refinement, based on maximum chord distance between mesh edge and geometric edge

*/
