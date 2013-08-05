/*
This section of the page provides templates for adding tables, and rows.

***************
* Empty Table *
***************

----------8<-----------------------------------------------------------8<----------
<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  0D Algorithms </b></big> </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  1D Algorithms </b></big> </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  2D Algorithms </b></big> </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  3D Algorithms </b></big> </td> </tr>

</tbody></table></center>
----------8<-----------------------------------------------------------8<----------

***************
* Table Entry *
***************

Color notes:
springgreen is for LGPL-compatible license - Bundled with Meshkit
skyblue is for LGPL-compatible license - External
goldenrod is for Propietary or GPL liscence External

Subpage note:
If an algorithm's description needs to be linked twice ONLY use \subpage once, for any occurance after the first, use \ref

----------8<-----------------------------------------------------------8<----------
<tr bgcolor=springgreen OR goldenrod OR skyblue>  <td> \endhtmlonly
<link to description i.e. \subpage desc_page>     \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \ref example_page>          \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \ref test_page>                \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>
----------8<-----------------------------------------------------------8<----------

*/

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









<br><br><big><big>  Mesh Generation Algorithms  </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  0D Algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_VertexMesher                   \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  1D Algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_EdgeMesher                     \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  2D Algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                            <td> \endhtmlonly
\subpage algorithm_CAMALTriAdvance                \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_CopyMesh                       \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_JaalQuad                       \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  3D Algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                            <td> \endhtmlonly
\subpage algorithm_CAMALTetMesher                 \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\ref algorithm_CopyMesh                           \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_ExtrudeMesh                    \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_EBMesh                         \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=skyblue>                              <td> \endhtmlonly
\subpage algorithm_NGTetMesher                    \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_OneToOneSwept                  \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_SCDMesh                        \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

</tbody></table></center>









<br><br><big><big>  Geometry Algorithms Algorithms  </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  2D Algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_AssyGen                        \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_CopyGeom                       \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  3D Algorithms </b></big> </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\ref algorithm_AssyGen                            \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\ref algorithm_CopyGeom                           \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_MeshOpTemplate                 \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

</tbody></table></center>









<br><br><big><big>  Mesh Modification Algorithms  </big></big><br><br>

<center><table width=80% bgcolor=black style="text-align: center;" border="0" cellpadding="5" cellspacing="2"><tbody>

<tr bgcolor=darkgray> <td>
Algorithm</td><td>
Related Examples</td><td>
Related Tests</td><td>
Brief Info</td></tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  2D Algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                            <td> \endhtmlonly
\subpage algorithm_MergeMesh                      \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=skyblue>                              <td> \endhtmlonly
\subpage algorithm_MesquiteOpt                    \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=springgreen>                          <td> \endhtmlonly
\subpage algorithm_QslimMesher                    \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=white> <td colspan="4"> <big><b>  3D Algorithms </b></big> </td> </tr>

<tr bgcolor=goldenrod>                            <td> \endhtmlonly
\ref algorithm_MergeMesh                          \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

<tr bgcolor=skyblue>                              <td> \endhtmlonly
\ref algorithm_MesquiteOpt                        \htmlonly </td> <td> \endhtmlonly
<link to example i.e. \\ref example_page>         \htmlonly </td> <td> \endhtmlonly
<link to test i.e. \\ref test_page>               \htmlonly </td> <td> \endhtmlonly
<a brief description, objective, short>           \htmlonly </td> </tr>

</tbody></table></center>

\endhtmlonly

*/
