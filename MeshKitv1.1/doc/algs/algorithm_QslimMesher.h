/*!
 \page algorithm_qslimmesher QslimMesher

<b>Name:</b> QslimMesher

<b>External dependencies:</b> (none) Adapted from http://mgarland.org/software/qslim10.html

<b>Input:</b> 2D ModelEnt's, already meshed with triangles

<b>Output:</b> Mesh vertices, triangles, in the original moab entity set

<b>Interface(s) used:</b> MOAB

<b>Setup:</b>

set_options() QslimOptions is used to set options specific to Qslim.

<b>Notes:</b>

  This decimation algorithm uses edge collapse as a primary
  simplification method. A cost for each possible edge collapse is
  established using quadric error concept (http://mgarland.org/research/quadrics.html)
*/
