/*!
 \page algorithm_VertexMesher VertexMesher

<b>Name:</b> VertexMesher

<b>External dependencies:</b> (none)

<b>Input:</b> 0D ModelEnt's

<b>Output:</b> Mesh vertices, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

None.

<b>Notes:</b>

This class implements vertex meshing.  A single instance of this class is created, since vertex meshing
is trivial.  Normally, applications need not interact with this mesher, since it gets created automatically
in the EdgeMesh setup() function.  The singleton instance can be accessed from the MeshKit instance by 
calling MeshKit::vertex_mesher().
*/
