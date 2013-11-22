/*!
 \page algorithm_extrudemesh ExtrudeMesh

<b>Name:</b> ExtrudeMesh

<b>External dependencies:</b> (none)

<b>Input:</b> 1D or 2D ModelEnt's, <em>already meshed</em>

<b>Output:</b> Extruded entities, and sometimes additional sets (see Notes below)

<b>Interface(s) used:</b> MOAB

<b>Setup:</b>

None.

<b>Notes:</b>

Extrudes mesh entities, along a straight or rotated path.

This algorithm accepts specification of copy, expand and extrude sets (see \ref copyexpandextrude_sets); 
see ExtrudeMesh header for details on specifying these sets.

*/
