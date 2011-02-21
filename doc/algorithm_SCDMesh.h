/*!
 \page algorithm_SCDMesh Structured Block Mesher

<b>Name:</b> SCDMesh

<b>External dependencies:</b> (none)

<b>Input:</b> 3D ModelEnt's

<b>Output:</b> Mesh vertices, hexes, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

None.

<b>Notes:</b>

This class generates a simple rectangular structured mesh, sized to completely surround the input 
ModelEnt(s).  The number of divisions in each parametric direction can be specified using member functions.
By default, the block is axis-aligned; however, this can also be changed by calling member functions.

*/
