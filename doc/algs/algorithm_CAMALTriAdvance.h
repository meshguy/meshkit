/*!
 \page algorithm_camaltriadvance CAMAL TriAdvance advancing front mesher

<b>Name:</b> CAMALTriAdvance

<b>External dependencies:</b> CAMAL Library (https://cubit.sandia.gov/public/camal.html)

<b>Input:</b> 2D ModelEnt's

<b>Output:</b> Mesh vertices, triangles, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

Calls MeshOp::setup_boundary().

<b>Notes:</b>

This class is the interface to the CAMAL TriAdvance algorithm.  This is a real-space, advancing front 
triangle meshing algorithm, implemented as part of the CUBIT mesh generation toolkit.  Currently only
isotropic size-based generation is supported.
*/
