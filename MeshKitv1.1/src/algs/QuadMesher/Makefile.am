# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libJaalQuadMesher.la
         
SUBDIRS = 

# The directory in which to install headers
libJaalQuadMesher_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libJaalQuadMesher_la_SOURCES = \
	basic_math.hpp \
	BinaryTree.cpp \
	BinaryTree.hpp \
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
	Mesh.cpp \
	Mesh.h \
	Mesh.hpp \
	mesh_io.cpp \
	MeshOpt.cpp \
	off_io.cpp \
	OneDefectRemeshing.cpp \
	QuadCleanUp.cpp \
	QuadCleanUp.hpp \
	QuadPatches.cpp \
        QuadMesh.cpp \
	Singlet.cpp \
	SwapEdges.hpp \
	SwapTriEdges.cpp \
	SwapQuadEdges.cpp \
	tfiblend.cpp \
	tfiblend.hpp \
	Tri2Quad.cpp \
	Tri2Quad.hpp \
	tri_io.cpp \
	circumcenter.cpp \
	circumcenter.hpp \
	GeomPredicates.cpp \
	GeomPredicates.hpp \
	MeshRefine2D.hpp \
	MeshRefine2D.cpp \
	MeshQuality.cpp  \
	StopWatch.hpp   \
	vtk_io.cpp      \
	xml_io.cpp     \
        cubit_io.cpp 
 
  

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
