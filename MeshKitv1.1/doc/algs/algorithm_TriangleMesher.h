/*!
 \page algorithm_trianglemesher TriangleMesher

<b>Name:</b> TriangleMesher

<b>External dependencies:</b> Triangle Library from http://www.cs.cmu.edu/~quake/triangle.html

<b>Input:</b> direction, input file

<b>Output:</b> 2D mesh

<b>Interface(s) used:</b> MOAB

<b>Setup:</b> 


<b>Notes:</b>
This TriangleMesher will mesh a set of nodes, and add just the triangles.
meshing is done in a plane controlled by a separate option (direction 1, 2, or 3)
default is 3, which means xy plane
options passed now are 'pc', which means that the convex hull is generated.
Triangles with one edge larger than a user suppled value (default HUGE) are removed
from the output.
Future improvements should allow for refinement of an existing mesh and constraining edges.
The user should keep in mind that everything happens in a plane.

MeshKit can link with static and shared-library for Triangle. Install Triangle using Ubuntu/OSX package manager.
If needed use commands:
"ar cru libtriangle.a triangle.o" "libtool --mode=compile gcc -Wall -c triangle.c" and/or libtool --mode=link gcc -rpath=/usr/local/lib -o libtriangle.la triangle.lo

*/
