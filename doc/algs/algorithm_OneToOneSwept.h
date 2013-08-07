/*!
 \page algorithm_onetooneswept One To One Sweeper

<b>Name:</b> OneToOneSwept

<b>External dependencies:</b> (none)

<b>Input:</b> 3D ModelEnt's

<b>Output:</b> Mesh vertices, quads, hexes, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

Source and linking surfaces checked for quad meshing algorithm, and for proper orientation of the mesh
to allow sweeping.

<b>Notes:</b>

This class sweep-meshes a volume, using single-source to single-target meshing.  Roca's algorithm is used
to compute interior vertex positions along the sweep.

<em>This algorithm sometimes generates inverted mesh for non-simply-connected source/target surfaces, especially
for large variations in the relative positions or sizes of internal loops between source and target surfaces.</em>

*/
