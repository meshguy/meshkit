/*!
  \page faq Frequently Asked Questions
 - \ref newdir "I'm implementing a new algorithm; should I put the source files directly in src/algs, or add a new subdirectory below that?"\n
 - \ref buildingnewdir "How do I add a new subdirectory, e.g. a new mesher with multiple files, to the &MeshKit source?"

\anchor newdir <b>I'm implementing a new algorithm; should I put the source files directly in src/algs, or add a new 
subdirectory below that?</b>\n
In general, if there are fewer than about 4 new files, put them in algs, otherwise put them in a new subdirectory.

\anchor buildingnewdir <b>How do I add a new subdirectory, e.g. a new mesher with multiple files, to the %MeshKit source?</b>\n
Let's say you have a new algorithm, myalg, that is implemented in several files.  If that number of files is
greater than about 4, we recommend you put it in its own subdirectory, under the src/algs directory in the %MeshKit
source.  To integrate this with the %MeshKit autotools-based build process, do the following steps:
 -# Create your subdirectory src/algs/myalg and its include directory src/algs/myalg/meshkit
 -# In src/algs/myalg/meshkit, put any include files you expect the application to see;
these are the headers that contain the API to your new mesher.
 -# Put the remaining files, both header and source files, in src/algs/myalg.
 -# In the top-level configure.ac file, just before the end of the file, add a line
to create the src/algs/myalg/Makefile (make it look like the others).
 -# In src/algs/Makefile.am, add myalg to the list of directories assigned to SUBDIRS
 -# In src/algs/Makefile.am, add your library to the LIBADD list for
     libMeshkitalgs.la
 -# Add your directory to the AM_CPPFLAGS in src/algs/Makfile.am so that
     you can include your algorithm in src/algs/register_algs.cpp
 -# Add an entry to src/algs/register_algs.cpp to register your algorithm
     with MKCore. (See \ref newmeshop).
 -# In src/algs/myalg, create a Makefile.am.  It is probably easiest to
    copy an existing one from elsewhere (e.g. some other subdir of algs/).
 -# Append to AM_CPPFLAGS the list of places that you expect the compiler
    to search for header files.  Keep in mind that for building in a separate
    directory from the source directory to work you will need to prefix
    paths with \$(srcdir)/.
 -# In the variable nobase_libMeshKitmyalg_la_include_HEADERS, substitute the names of all the files you put in
src/algs/myalg/meshkit
 -# Put the names of all the files you put in src/algs/myalg in libMeshKitmyalg_la_SOURCES (including
the header files there)
 -# In the top-level %MeshKit directory, run autoreconf -fi
 -# Run make; the code in myalg should build

Top: \ref index Prev: \ref index Next: \ref index
  
 */
