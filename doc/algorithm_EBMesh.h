/*!
 \page algorithm_EBMesh Embedded Boundary Mesher

<b>Name:</b> EBMesh

<b>External dependencies:</b> (none)

<b>Input:</b> 3D ModelEnt's, and a rectangular region of structured mesh (e.g. created by SCDMesh)

<b>Output:</b> Tags set on hexes in structured mesh, indicating inside, boundary, and outside elements

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

Looks on the ModelEnt for an assigned SCDMesh operation, and if none is found, assigns one.

<b>Notes:</b>

This class implements embedded boundary mesh generation, based on ray tracing implemented in MOAB's 
Hierarchical OBB tree.

*/
