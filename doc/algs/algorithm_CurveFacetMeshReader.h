/*!
 \page algorithm_curvefacetmeshreader CurveFacetMeshReader

<b>Name:</b> CurveFacetMeshReader

<b>External dependencies:</b> (none)

<b>Input:</b> 1D ModelEnt's

<b>Output:</b> Mesh edges & vertices, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

None.

<b>Notes:</b>


This class uses the facets generated for the visualization of solid models in typical CAD 
software to represent a geometric curve as mesh. Upon execution, this class will call for
the facet data and store it as part of the ModelEnt's mesh.

*/
