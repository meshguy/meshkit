/*!
\page meshkitalgorithms MeshKit Algorithms
Any algorithm in section "Undocumented Algorithms" eventually need place in appropiate other tables.\n
Any cell that has "(cell needs review)" needs it's content and the content on the other side of links reviewed.

\htmlonly

<center>
<table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody><tr>
<td bgcolor=darkgray> Key: </td>
<td bgcolor=springgreen> LGPL-compatible license <br> Bundled with MeshKit </td>
<td bgcolor=skyblue> LGPL-compatible license <br> External </td>
<td bgcolor=goldenrod> Proprietary or GPL license <br> External </td>
</tr></tbody></table></center>


<br><br><big><big> Mesh Generation Algorithms </big></big><br><br>

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
meshes vertices                                                                               \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 1D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_edgemesher "EdgeMesher"                                                    \htmlonly </td> <td> \endhtmlonly
\ref example_edgemesher.cpp                                                                   \htmlonly </td> <td> \endhtmlonly
\ref test_edgemesher.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
meshes curves                                                                                 \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 2D algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camaltriadvance "CAMALTriAdvance"                                          \htmlonly </td> <td> \endhtmlonly
\ref example_camaltriadvance.cpp                                    \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camaltriadvance.cpp                                                                 \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_copymesh "CopyMesh"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_copymesh.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_copymesh.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_jaalquad "JaalQuad"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_jaalquad.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_jaalquad.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 3D algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camaltetmesher "CAMALTetMesher"                                            \htmlonly </td> <td> \endhtmlonly
\ref example_camaltetmesher.cpp                                     \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camaltetmesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\ref algorithm_copymesh "CopyMesh"                                                            \htmlonly </td> <td> \endhtmlonly
\ref example_copymesh.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_copymesh.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_extrudemesh "ExtrudeMesh"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_extrudemesh.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_extrudemesh.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_ebmesh "EBMesh"                                                            \htmlonly </td> <td> \endhtmlonly
\ref example_ebmesh.cpp                                             \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_ebmesh.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=skyblue>                                                                          <td> \endhtmlonly
\subpage algorithm_ngtetmesher "NGTetMesher"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_ngtetmesher.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_ngtetmesher.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_onetooneswept "OneToOneSwept"                                              \htmlonly </td> <td> \endhtmlonly
\ref example_onetooneswept.cpp                                                                \htmlonly </td> <td> \endhtmlonly
\ref test_onetooneswept.cpp                                                                   \htmlonly </td> <td> \endhtmlonly
sweeper                                                                                       \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_scdmesh "SCDMesh"                                                          \htmlonly </td> <td> \endhtmlonly
\ref example_scdmesh.cpp                                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_scdmesh.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

</tbody></table></center>


<br><br><big><big> Geometry Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 2D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_assygen "AssyGen"                                                          \htmlonly </td> <td> \endhtmlonly
\ref example_assygen.cpp                                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_assygen.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_copygeom "CopyGeom"                                                        \htmlonly </td> <td> \endhtmlonly
\ref example_copygeom.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_copygeom.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 3D algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\ref algorithm_assygen "AssyGen"                                                              \htmlonly </td> <td> \endhtmlonly
\ref example_assygen.cpp                                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_assygen.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\ref algorithm_copygeom "CopyGeom"                                                            \htmlonly </td> <td> \endhtmlonly
\ref example_copygeom.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_copygeom.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_meshoptemplate "MeshOpTemplate"                                            \htmlonly </td> <td> \endhtmlonly
\ref example_meshoptemplate.cpp                                     \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_meshoptemplate.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

</tbody></table></center>


<br><br><big><big> Mesh Modification Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 2D algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_mergemesh "MergeMesh"                                                      \htmlonly </td> <td> \endhtmlonly
(none)                                                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
(none)                                                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=skyblue>                                                                          <td> \endhtmlonly
\subpage algorithm_mesquiteopt "MesquiteOpt"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_mesquiteopt.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_mesquiteopt.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_qslimmesher "QSlimMesher"                                                  \htmlonly </td> <td> \endhtmlonly
\ref example_qslimmesher.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_qslimmesher.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b> 3D algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\ref algorithm_mergemesh "MergeMesh"                                                          \htmlonly </td> <td> \endhtmlonly
(none)                                                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
(none)                                                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=skyblue>                                                                          <td> \endhtmlonly
\ref algorithm_mesquiteopt "MesquiteOpt"                                                      \htmlonly </td> <td> \endhtmlonly
\ref example_mesquiteopt.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_mesquiteopt.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

</tbody></table></center>

<br><br><big><big> Undocumented Algorithms </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b> XD Algorithms </b></big> </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_fbgeom "FBGeom"                                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_fbgeom.cpp                                             \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_fbgeom.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_intassign "intassign"                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_intassign.cpp                                          \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_intassign.cpp                                                                       \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camalmb2 "CAMALMB2"                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_camalmb2.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camalmb2.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camalmbgeom "CAMALMBGeom"                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_camalmbgeom.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camalmbgeom.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_mbgeomop "MBGeomOp"                              \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_mbgeomop.cpp                                           \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_mbgeomop.cpp                                                                        \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camalmbgeomvar "CAMALMBGeomVar"                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_camalmbgeomvar.cpp                                     \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camalmbgeomvar.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_mbsplitop "MBSplitOp"                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_mbsplitop.cpp                                          \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_mbsplitop.cpp                                                                       \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_mbvolop "MBVolOp"                                \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_mbvolop.cpp                                            \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_mbvolop.cpp                                                                         \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_parallelmesher "ParallelMesher"                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_parallelmesher.cpp                                     \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_parallelmesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=goldenrod>                                                                        <td> \endhtmlonly
\subpage algorithm_camalpaver "CamalPaver"                          \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_camalpaver.cpp                                         \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_camalpaver.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
Unstructured Quad Mesher                                            \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_postbl "PostBL"                                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_postbl.cpp                                             \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_postbl.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_setpnt2quad "SetPnt2Quad"                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_setpnt2quad.cpp                                        \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_setpnt2quad.cpp                                                                     \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=springgreen>                                                                      <td> \endhtmlonly
\subpage algorithm_tfimapping "TFIMapping"                          \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_tfimapping.cpp                                                                   \htmlonly </td> <td> \endhtmlonly
\ref test_tfimapping.cpp                                                                      \htmlonly </td> <td> \endhtmlonly
structured quad mesher                                                                        \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_trianglemesher "TriangleMesher"                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_trianglemesher.cpp                                     \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_trianglemesher.cpp                                                                  \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tr bgcolor=BBFFBB>                                                                           <td> \endhtmlonly
\subpage algorithm_volice "volIce"                                  \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref example_volice.cpp                                             \n(cell needs review)     \htmlonly </td> <td> \endhtmlonly
\ref test_volice.cpp                                                                          \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)                          \n(cell needs review)     \htmlonly </td> </tr>

<tbody></table></center>

\endhtmlonly


*/
