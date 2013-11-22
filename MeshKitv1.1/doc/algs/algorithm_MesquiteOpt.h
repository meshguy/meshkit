/*!\page algorithm_mesquiteopt MesquiteOpt

<b>Name:</b> MesquiteOpt

<b>External dependencies:</b> Mesquite (http://www.cs.sandia.gov/optimization/knupp/Mesquite.html)

<b>Input:</b> 2D or 3D ModelEnts, <em>already meshed</em>

<b>Output:</b> The input mesh with vertex positions adjusted in improve mesh quality

<b>Interface(s) used:</b> MOAB, iMesh, iGeom, Mesquite

<b>Setup:</b>

None required.  Default is to run ShapeImprovementWraper (or ShapeImprover 
depending on the version of Mesquite) with the mesh boundary fixed.

\c smooth_with_free_boundary() can be called to change to a "free" smooth, where
vertices on the boundary of the mesh may also be moved.  This option requires
child geomery defining the boundary (with the obvious exception of smoothing a 
surface that is a topological sphere) so that the optimization problem is
properly constrained.

\c set_mesh_op() may be called to specify an alogithm other than the default
shape improver to run.

*/

