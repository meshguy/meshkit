/*!
 \page algorithm_camaltetmesher CAMAL Tet mesher

<b>Name:</b> CAMALTetMesher

<b>External dependencies:</b> CAMAL Library (https://cubit.sandia.gov/public/camal.html)

<b>Input:</b> 3D ModelEnt's

<b>Output:</b> Mesh vertices, tetrahedra, committed to ModelEnt

<b>Interface(s) used:</b> MOAB

<b>Setup:</b>

Looks on the ModelEnt for an assigned tri mesher, and if none is found, assigns one from the default list
of tri meshing algorithms registered with MeshKit.

<b>Notes:</b>

This class is the interface to the CAMAL Tet meshing algorithm, which is itself a wrapper for the
Simulog tetrahedral meshing algorithm.  Currently only isotropic size-based generation is supported.
*/
