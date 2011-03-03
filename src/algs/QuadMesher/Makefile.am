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
	basic_math.hpp \
	BinaryTree.cpp \
	BinaryTree.hpp \
	Bridges.cpp \
	cycle.hpp \
	DijkstraShortestPath.cpp \
	DijkstraShortestPath.hpp \
	Doublet.cpp \
	DualGraph.cpp \
	DualGraph.hpp \
	FastArea.hpp \
	FaceClose.cpp \
	JaalMoabConverter.cpp \
	JaalMoabConverter.hpp \
	laplace.cpp \
	matching.cpp \
	Mesh.cpp \
	Mesh.h \
	Mesh.hpp \
	meshclean.cpp \
	mesh_io.cpp \
	MeshOpt.cpp \
	MeshRefine2D.hpp \
	meshViewer.hpp \
	ObjectPool.hpp \
	off_io.cpp \
	OneDefectRemeshing.cpp \
	QuadCleanUp.cpp \
	QuadCleanUp.hpp \
	QuadPatches.cpp \
        QuadMesh.cpp \
	Singlet.cpp \
	SwapQuadEdge.cpp \
	SwapTriEdge.cpp \
	SwapTriEdge.hpp \
	tfiblend.cpp \
	tfiblend.hpp \
	Tri2Quad.cpp \
	Tri2Quad.hpp \
	tri_io.cpp \
	Tunnels.cpp \
	circumcenter.cpp \
	circumcenter.hpp \
	geomPredicates.cpp \
	geomPredicates.hpp \
	MeshRefine2D.cpp \
	QuadQuality.cpp  \
	vtk_io.cpp
  

# The list of header files which are to be installed
nobase_libJaalQuadMesher_la_include_HEADERS = \
        meshkit/QuadMesh.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../../core \
              -I$(srcdir)/../../lemon \
              -I$(srcdir)/../../utils \
              -I$(top_builddir) \
              -I$(top_srcdir)/utils 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
