/*!
\page meshkitalgorithms MeshKit Algorithms
\htmlonly

<center>
<table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody><tr>
<td bgcolor=darkgray> Key: </td>
<td bgcolor=springgreen> LGPL-compatible license <br> Bundled with MeshKit </td>
<td bgcolor=skyblue> LGPL-compatible license <br> External </td>
<td bgcolor=goldenrod> Proprietary or GPL license <br> External </td>
</tr></tbody></table></center>


<br><br><big><big> Unstructured Mesh Generation Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 0D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_vertexmesher "VertexMesher"                                                \htmlonly </td> <td> \endhtmlonly
(none)                                                                                        \htmlonly </td> <td> \endhtmlonly
(none)                                                                                        \htmlonly </td> <td> \endhtmlonly
Meshes vertices                                                                               \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 1D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_edgemesher "EdgeMesher"                                                    \htmlonly </td> <td> \endhtmlonly
\ref example_edgemesher.cpp                                                                   \htmlonly </td> <td> \endhtmlonly
\ref test_edgemesher.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
Meshes curves                                                                                 \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 2D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_tfimapping "TFIMapping"                          \n     \htmlonly </td> <td> \endhtmlonly
\ref example_tfimapping.cpp example_basic.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
\ref test_tfimapping.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
Mapping algorithm for quads                                                                        \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_jaalquad "JaalQuad"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_jaalquad.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_jaalquad.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Quad mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camaltriadvance "CAMALTriAdvance"                                          \htmlonly </td> <td> \endhtmlonly
\ref example_camaltriadvance.cpp                                    \n     \htmlonly </td> <td> \endhtmlonly
\ref test_camaltriadvance.cpp                                                                 \htmlonly </td> <td> \endhtmlonly
Advancing-front tri mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camalpaver "CamalPaver"                          \n     \htmlonly </td> <td> \endhtmlonly
\ref example_camalpaver.cpp                                         \n     \htmlonly </td> <td> \endhtmlonly
\ref test_camalpaver.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
Unstructured quad mesher                                            \n     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                           <td> \endhtmlonly
\subpage algorithm_trianglemesher "TriangleMesher"                  \n     \htmlonly </td> <td> \endhtmlonly
\ref example_trianglemesher.cpp                                     \n     \htmlonly </td> <td> \endhtmlonly
\ref test_trianglemesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
Delaunay-based tri mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 3D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_onetooneswept "OneToOneSwept"                                              \htmlonly </td> <td> \endhtmlonly
\ref example_onetooneswept.cpp example_basic.cpp                                                               \htmlonly </td> <td> \endhtmlonly
\ref test_onetooneswept.cpp                                                                   \htmlonly </td> <td> \endhtmlonly
Sweeper creates an all hex mesh.                                                                        \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camaltetmesher "CAMALTetMesher"                                            \htmlonly </td> <td> \endhtmlonly
\ref example_camaltetmesher.cpp                                     \n     \htmlonly </td> <td> \endhtmlonly
\ref test_camaltetmesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
Tet mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=skyblue>                                                                          <td> \endhtmlonly
\subpage algorithm_ngtetmesher "NGTetMesher"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_ngtetmesher.cpp                                        \n     \htmlonly </td> <td> \endhtmlonly
\ref test_ngtetmesher.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
Tet mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_parallelmesher "ParallelMesher"                  \n     \htmlonly </td> <td> \endhtmlonly
\ref example_parallelmesher.cpp                                     \n     \htmlonly </td> <td> \endhtmlonly
\ref test_parallelmesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
CAD-based parallel mesher                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_extrudemesh "ExtrudeMesh"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_extrudemesh.cpp                                        \n     \htmlonly </td> <td> \endhtmlonly
\ref test_extrudemesh.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
Extrudes 1D or 2D mesh.                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_meshoptemplate "MeshOpTemplate"                                            \htmlonly </td> <td> \endhtmlonly
\ref example_meshoptemplate.cpp example_basic.cpp                                    \n     \htmlonly </td> <td> \endhtmlonly
\ref test_meshoptemplate.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
A dummy algorithm template that creates geometry                         \n     \htmlonly </td> </tr>\

</tbody></table></center>



</tbody></table></center>


<br><br><big><big> Structured Mesh Generation Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_scdmesh "SCDMesh"                                                          \htmlonly </td> <td> \endhtmlonly
\ref example_scdmesh.cpp                                            \n     \htmlonly </td> <td> \endhtmlonly
\ref test_scdmesh.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
Simple rectangular structured mesh                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_ebmesh "EBMesh"                                                            \htmlonly </td> <td> \endhtmlonly
\ref example_ebmesh.cpp                                             \n     \htmlonly </td> <td> \endhtmlonly
\ref test_ebmesh.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
Embedded boundary mesher                          \n     \htmlonly </td> </tr>

<tbody></table></center>

<br><br><big><big> Mesh Modification Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=skyblue>                                                                          <td> \endhtmlonly
\subpage algorithm_mesquiteopt "MesquiteOpt"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_mesquiteopt.cpp                                        \n     \htmlonly </td> <td> \endhtmlonly
\ref test_mesquiteopt.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
Mesh quality optimizer                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_postbl "PostBL"                                  \n     \htmlonly </td> <td> \endhtmlonly
\ref example_postbl.cpp                                             \n     \htmlonly </td> <td> \endhtmlonly
\ref test_postbl.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
Add boundary layers to existing mesh                         \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_qslimmesher "QSlimMesher"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_qslimmesher.cpp                                        \n     \htmlonly </td> <td> \endhtmlonly
\ref test_qslimmesher.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
Mesh decimatator                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_intassign "IntervalAssignment"                            \n     \htmlonly </td> <td> \endhtmlonly
\ref example_intassign.cpp                                          \n     \htmlonly </td> <td> \endhtmlonly
\ref test_intassign.cpp                                                                       \htmlonly </td> <td> \endhtmlonly
Interval assignment algorithm                         \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_copymesh "CopyMesh"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_copymesh.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_copymesh.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Copy moves an existing mesh                         \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                        <td> \endhtmlonly
\subpage algorithm_mergemesh "MergeMesh"                                                      \htmlonly </td> <td> \endhtmlonly
\ref example_copymesh.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_copymesh.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Merges duplicate nodes                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_copygeom "CopyGeom"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_copygeom.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_copygeom.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Copy moves existing geometry                          \n     \htmlonly </td> </tr>

</tbody></table></center>

<br><br><big><big> Nuclear Reactor Modeling Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_assygen "AssyGen"                                                          \htmlonly </td> <td> \endhtmlonly
\ref example_assygen.cpp                                            \n     \htmlonly </td> <td> \endhtmlonly
\ref test_assygen.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
Creates nuclear reactor assemblies                          \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_coregen "CoreGen"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_coregen.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_coregen.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Creates nuclear reactor core from assembly meshes                         \n     \htmlonly </td> </tr>
</tbody></table></center>

<br><br><big><big> Mesh-based Geometry Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_mbgeomop "MBGeomOp"                              \n     \htmlonly </td> <td> \endhtmlonly
\ref example_mbgeomop.cpp                                           \n     \htmlonly </td> <td> \endhtmlonly
\ref test_mbgeomop.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
Creates geometrized surface mesh from manifold surface mesh                         \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_mbsplitop "MBSplitOp"                            \n     \htmlonly </td> <td> \endhtmlonly
\ref example_mbsplitop.cpp                                          \n     \htmlonly </td> <td> \endhtmlonly
\ref test_mbsplitop.cpp                                                                       \htmlonly </td> <td> \endhtmlonly
Edits - crops/split surface                         \n     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                           <td> \endhtmlonly
\subpage algorithm_mbvolop "MBVolOp"                                \n     \htmlonly </td> <td> \endhtmlonly
\ref example_mbvolop.cpp                                            \n     \htmlonly </td> <td> \endhtmlonly
\ref test_mbvolop.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
Creates a meshable BREP volume of interest                          \n     \htmlonly </td> </tr>
</tbody></table></center>



\endhtmlonly


*/
