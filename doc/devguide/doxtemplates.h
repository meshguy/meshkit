/** \page doxtemplates Doxygen Algorithm Templates

\tableofcontents

\section mkalgs MeshKit Algorithms Page

In all templates:
- Parentheses sections should be replaced
- Parentheses should be removed


\subsection mkalgs_table Table

If there is not a table where a newly created algorithm fits a new table can be inserted.

Below are some notes on using the template that follows after:
- A generic title for table should be made so other algorithms will fit
- Dimensions that are not represented in algorithm should be removed


\verbatim
<br><br><big><big>  (name of table)  </big></big><br><br>

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
\endverbatim


\subsection mkalgs_entry Entry

Each algorithm should have at least 1 entry on the Meshkit Algorithms page.

Below are some notes on using the template that follows after:
- Color must be changed.
  - "springreen" is for LGPL-compatible license; Bundled with MeshKit 
  - "skyblue" is for LGPL-compatible license; External 
  - "goldenrod" is for Proprietary or GPL license; External 
- No algorithm should have two "\subpage" links, use "\ref" for any past the first
- If multiple tests/examples exist:
  - Tests can be added on same line, and will wrap when neccesary
  - Newlines should only be added by adding "\n"

Template:
\verbatim
<tr bgcolor=(color)>                                  <td> \endhtmlonly
(link to description i.e. \subpage alg_page "page")   \htmlonly </td> <td> \endhtmlonly
(link to example i.e. \ref example_page)              \htmlonly </td> <td> \endhtmlonly
(link to test i.e. \ref test_page)                    \htmlonly </td> <td> \endhtmlonly
(a brief description; objective and short)            \htmlonly </td> </tr>
\endverbatim


\section descriptions Algorithm Description

For each algorithm that is created, an algorithm description page should be created.
The purpose of this page is to give users a brief look at the algorithm and it's requirements.

\verbatim
/*!
\page algorithm_(lowercase meshop name) (name of this algorithm)

<b>Name:</b> (name of this algorithm)

<b>External dependencies:</b> (none)

<b>Input:</b> (valid input to this algorithm)

<b>Output:</b> (output generated by this algorithm)

<b>Interface(s) used:</b> (selection of iGeom, iMesh, iRel, CGM, MOAB, LASSO)

<b>Setup:</b>

(any noteable calls during setup, or requirements for successful setup/execute belong here [Pre-Setup])

<b>Notes:</b>

(any notes on using the algorithm belong here)
*/
\endverbatim


\section tests Tests
At the top of every test there should be a doxygen section with some brief infromation.

Template:
\verbatim
/**
\file test_(lowercase meshop or feature).cpp
\test This file contains (some number of tests, i.e. "3") test(s). Algorithms tested include: (List of MeshOp's or features tested).
*/
\endverbatim


\section examples Examples

At the top of every example there should be a doxygen section with some brief information.
Examples should also contain many C-style comments in the code so new users can follow along.

Below are some notes on using the template that follows after:
- Sections with "(myalg)" should be replaced by
  - A lowercase name of MeshOp if appropiate
  - A lowercase short name as a description of example otherwise
- \\bug and \\warning can be removed if there are none
- \\image
  - If there is no input or output, the line should be replaced with "(none)"
  - You must add create and add the images to trunk/doc/figures/examples/
  - There is a shell script that will resize images (trunk/dox/resize.sh)
    - It has a default size built in appropiate for doxygen
    - Modifies the original image
    - Depends on the programs "convert", "expr", and "identify" being available


Template:
\verbatim
/**
\example example_(myalg).cpp

\section example_(myalg)_cpp_title <name of My Algorithm>

(some text explaining objectives of example)

\subsection example_(myalg)_cpp_in Input
\image html (myalg).in.jpg "(description of image)"

\subsection example_myalg_cpp_out Output
\image html (myalg).out.jpg "(description of image)"

\subsection example_(myalg)_cpp_inf Misc. Information
\author (author's name)
\date (date-of-creation)
\bug (message if there are any bugs with your example)
\warning (message if users should be warned about your example)

\subsection example_myalg_cpp_src Source Code
*/
\endverbatim

*/