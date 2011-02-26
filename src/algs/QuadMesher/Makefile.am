# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libJaalQuadMesher.la
         
SUBDIRS = 

# The directory in which to install headers
libJaalQuadMesher_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libJaalQuadMesher_la_SOURCES = \
	AdfrontClean.cpp  \
	BinaryTree.cpp \
	Bridges.cpp \
	DijkstraShortestPath.cpp \
	Doublet.cpp \
	DualGraph.cpp \
	FaceClose.cpp \
	laplace.cpp \
	matching.cpp \
	Mesh.cpp \
	meshclean.cpp \
	mesh_io.cpp \
	MeshOpt.cpp \
	off_io.cpp \
	OneDefectRemeshing.cpp \
	QuadCleanUp.cpp \
	QuadPatches.cpp \
	Singlet.cpp \
	SwapQuadEdge.cpp \
	SwapTriEdge.cpp \
	tfiblend.cpp \
	Tri2Quad.cpp \
	tri_io.cpp \
	Tunnels.cpp \
	circumcenter.cpp \
	geomPredicates.cpp \
	MeshRefine2D.cpp \
	QuadQuality.cpp  \
	vtk_io.cpp
  

# The list of header files which are to be installed
nobase_libJaalQuadMesher_la_include_HEADERS = \
	meshkit/Mesh.hpp   \
	meshkit/MeshRefine2D.hpp \
        meshkit/QuadMesher.hpp \
	meshkit/JaalMoabConverter.hpp \
	meshkit/QuadCleanUp.hpp  \
	meshkit/SwapTriEdge.hpp  \
	meshkit/Tri2Quad.hpp	

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core \
              -I$(srcdir)/../lemon \
              -I$(srcdir)/../utils \
              -I$(top_builddir) \
              -I$(top_srcdir)/utils 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
