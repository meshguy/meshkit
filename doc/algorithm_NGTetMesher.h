/*!
 \page algorithm_NGTetMesher Netgen Tet mesher

<b>Name:</b> NGTetMesher

<b>External dependencies:</b> \ref extern_netgen

<b>Input:</b> 3D ModelEnt's

<b>Output:</b> Mesh vertices, tetrahedra, committed to ModelEnt

<b>Interface(s) used:</b> MOAB

<b>Setup:</b>

Looks on the ModelEnt for an assigned tri mesher, and if none is found, assigns one from the default list
of tri meshing algorithms registered with MeshKit.

<b>Notes:</b>

This class is the interface to the Netgen Tet meshing algorithm.  Currently only isotropic size-based 
generation is supported.
*/
