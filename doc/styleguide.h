/*!
  \page styleguide Coding Style Guide
Code developed in %MeshKit should follow the coding styles described here.  Any deviations from this style
guide will result in severe berating and other verbal abuse.

\section dirstructure MeshKit Directory Structure
%MeshKit source code is organized in the following directory structure: \n
 - data: This directory contains geometry and mesh data files used in unit tests.
 - doc: Documentation is put here, along with the input file for Doxygen.  Most %MeshKit documentation is doxygen-processed.
 - examples: Examples of %MeshKit usage, both small and large.  These programs are not meant to be used as unit tests, but
     rather as further documentation on %MeshKit usage.
 - src
   - core: %MeshKit core classes; these classes should not change very often compared to the rest of %MeshKit.
   - algs: Meshing algorithms, both those that generate as well as modify mesh.
     - <myalg_subdir>: Directories below src/algs should only be used in cases where an algorithm has more than a few
         (roughly speaking, > 2) source files.
   - utils: Utility classes, used in algorithms and elsewhere.  These classes should mostly be lightweight and therefore
            inexpensive to construct and destruct.
 - test: All unit test programs should go below this directory.
   - algs: Unit tests for the meshing algorithms in %MeshKit.
   - core: Unit tests for the core classes in %MeshKit.
If you're designing a new class or other code for %MeshKit and are not sure where to put it, try to find something similar
and put it there.  Otherwise, email the %MeshKit email list for pointers.  <em> In general, you should not need to create new
subdirectories in the %MeshKit source code, except when implementing a new algorithm with more than about 2 files.</em>

\section sourcestyle Source Code Style
%MeshKit code should abide by the following general rules:
 - Names:
   - Class names should be in the CamelBack style, e.g. EdgeMesh or VertexMesher.
   - Class member variables should be camelBack, e.g. EdgeMesh::schemeType; each member variable, e.g. int memberVariable, 
   should have set/get functions void member_variable(int newval) and int member_variable(), respectively.
   - Enumeration values should be all captitalized, with underscores avoided if possible (the enumeration name indicates
     the general purpose of the enumeration, so e.g. we use EQUAL, not EQUAL_MESH)
 - Source code should not contain tabs or MS-DOS newlines; tabs and other indentations should be set to a width of 2 spaces.
   For general tips on how to set your editor for this, see the %MeshKit-dev discussion starting with <a href="https://lists.mcs.anl.gov/mailman/private/meshkit-dev/2011/000519.html">this message</a>.
 - Each class header should be fully commented; that includes:
   - A \\file comment block at the top of the file; DO NOT include things like Author and Date blocks; this stuff is available
     from subversion if we really need to know.
   - A \\class comment block, formatted like those in the %MeshKit core classes.  THE FIELDS AFTER THE CLASS NAME ARE VERY IMPORTANT,
     as they tell developers how to include the class declaration in their code.  This information goes into the "Detailed
     Description" portion of the class documentation.  This block should include any features or limitations of the class.
     Eventually, at least for meshing classes, we'll impose some standard boilerplate that each meshing class should use.
     Until then, please keep this block to around a paragraph.
   - Each function in both the public and private interfaces should be commented, INCLUDING ANY ARGUMENTS AND RETURN VALUES.
     See the %MeshKit core classes for examples of how to format these comments.  As a rule of thumb, your code should run through
     Doxygen without generating any warnings; in fact, Doxygen is sometimes helpful at pointing out inconsistencies in your
     class declaration.
 - Developers should avoid using #include in header files, as they propagate dependencies more widely than necessary.  The only
   cases where other includes are needed are to import the declaration for a parent class, and to declare types used as
   non-pointer and non-reference function arguments.  In most cases, a forward-declaration statement (e.g. 'class MKCore') 
   will suffice.
   

\section commits Making Repository Commits
As a general rule, developers should update frequently, and commit changes often.  However, the repository should always remain
in a state where the code can be compiled.  Most of the time, the code should also successfully execute "make check" run from the
top-level directory.  If you commit code that violates this principal, it should be your first priority to return the repository
code to a compilable state, and your second priority to make sure "make check" runs without errors.

Commits to the repository should also come with a non-trivial, useful, non-verbose log message.  Oftentimes the best way to generate
this message is to run 'svn diff > diffs', and edit the diffs file to remove specific line changes but include a comment on 
each file that changed.  Many times it is helpful to state that 'make check runs successfully' at the end of the log message.
Although it would be possible and many software projects do it, we prefer not to force successful execution of the test suite 
before every commit.  Developers should make every effort to avoid having to impose this constraint, by running a make check
before every commit.

Top: \ref index Prev: \ref devinterfaces Next: \ref faq
  
 */
