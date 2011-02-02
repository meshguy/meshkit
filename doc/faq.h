/*!
  \page faq Frequently Asked Questions
 - \ref newdir "I'm implementing a new algorithm; should I put the source files directly in src/algs, or add a new subdirectory below that?"\n
 - \ref buildingnewdir "How do I add a new subdirectory, e.g. a new mesher with multiple files, to the %%MeshKit source?"

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
 -# In src/algs/myalg, create a Makefile.am, by copying the one from src/algs/Makefile.am
 -# Everywhere you see lib%MeshKitalgs, substitute your own library name, e.g. lib%MeshKitmyalg
 -# In the variable nobase_lib%MeshKitalgs_la_include_HEADERS, substitute the names of all the files you put in
src/algs/myalg/meshkit
 -# Put the names of all the files you put in src/algs/myalg in nobase_lib%MeshKitalgs_la_include_SOURCES (including
the header files there)
 -# Add the new directory src/algs/myalg to the AM_CPPFLAGS in src/algs/myalg/Makefile.am (with a -I in front of it)
 -# In the top-level %MeshKit directory, run autoreconf -fi
 -# Run make; the code in myalg should build

Top: \ref index Prev: \ref index Next: \ref index
  
 */
