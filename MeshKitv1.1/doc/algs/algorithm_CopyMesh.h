/*!
 \page algorithm_copymesh CopyMesh

<b>Name:</b> CopyMesh

<b>External dependencies:</b> (none)

<b>Input:</b> 1D, 2D, or 3D ModelEnt's, <em>already meshed</em>

<b>Output:</b> Copied entities, and sometimes additional sets (see Notes below)

<b>Interface(s) used:</b> MOAB

<b>Setup:</b>

None.

<b>Notes:</b>

Copies mesh entities.  Copied entities can optionally retain a tag pointing to the resulting copy.

This algorithm accepts specification of copy and expand sets (see \ref copyexpandextrude_sets); 
see CopyMesh header for details on specifying these sets.

*/
