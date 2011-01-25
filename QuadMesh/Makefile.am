# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
lib_LTLIBRARIES = libQuadMesh.la

# The directory in which to install headers
libMKalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libQuadMesh_la_SOURCES = \
	BinaryTree.cc \
	DualGraph.cc  \
	Mesh.cc       \
	QuadCleanUp.cc \
	Tri2Quad.cc \
	laplace.cc \
	mesh_io.cc \
	off_io.cc


# The list of header files which are to be installed
libQuadMesh_la_include_HEADERS = \
	BinaryTree.h  \
	DualGraph.h  \
	Mesh.h  \
	QuadCleanUp.h  \
	SimpleArray.h  \
	Tri2Quad.h

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)

if HAVE_MOAB
CXXFLAGS += -DMOAB
endif

