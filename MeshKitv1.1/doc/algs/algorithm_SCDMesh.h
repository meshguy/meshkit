/*!
 \page algorithm_scdmesh Structured Block Mesher

<b>Name:</b> SCDMesh

<b>External dependencies:</b> (none)

<b>Input:</b> 3D ModelEnt's

<b>Output:</b> Mesh vertices, hexes, committed to ModelEnt

<b>Interface(s) used:</b> iGeom, MOAB

<b>Setup:</b>

None.

<b>Notes:</b>

This class generates a simple rectangular structured mesh, sized to completely 
surround the input ModelEnt(s). Three options are defined for structured grid
generation that are defined below: mesh representation, grid size definition,
axis type.

<b>Mesh Representation Options:</b>

By default, a full representation of the mesh is created where every vertex is 
defined and every hexahedral element constructed from connectivity arrays. 
Alternatively, MOAB's structured mesh interface can be used to create the mesh,
resulting in a more compact representation and smaller memory usage. These 
options can be accessed through the C++ interface with the set_interface_scheme()
member function of the SCDMesh mesh operation as follows:

SCDMeshInstance->set_interface_scheme(SCDMesh::full);

or

SCDMeshInstance->set_interface_scheme(SCDMesh::scd);


<b>Grid Sizing Options:</b>

By default, grid sizing is defined by setting the number of coarse divisions
in the ijk axes and the number of fine divisions in each of those coarse
divisions. Optionally, a full array for each ijk axis can be defined to 
provided any arbitrary non-uniform spacing. These options can be accessed 
through the C++ interface with the set_grid_scheme() member function of the 
SCDMesh mesh operation as follows:

SCDMeshInstance->set_grid_scheme(SCDMesh::cfMesh);
int coarse_i_divisions;
std::vector<int> fine_i_divisions;
SCDMeshInstance->set_coarse_i_grid(coarse_i_divisions);
SCDMeshInstance->set_fine_i_grid(fine_i_divisions);

or

SCDMeshInstance->set_grid_scheme(SCDMesh::vtxMesh);
std::vector<double> i_coords;
SCDMeshInstance->set_i_coordinates(i_coords);

<b>Axis Type Options:</b>

By default, the grid is aligned on the xyz Cartesian axes. Optionally, the grid
can be oriented to the primary geometric axes of the model entity being meshed.
These options can be accessed through the C++ interface with the set_axis_scheme()
member function of the SCDMesh mesh operation as follows:

SCDMeshInstance->set_axis_scheme(SCDMesh::cartesian);

or 

SCDMeshInstance->set_axis_scheme(oriented);
*/
