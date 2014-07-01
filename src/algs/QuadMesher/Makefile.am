# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libJaalQuadMesher.la
         
SUBDIRS = 

# The directory in which to install headers
libJaalQuadMesher_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libJaalQuadMesher_la_SOURCES = \
	BinaryTree.cpp \
	DijkstraShortestPath.cpp \
	Doublet.cpp \
	DualGraph.cpp \
	FaceClose.cpp \
	JaalMoabConverter.cpp \
	laplace.cpp \
	Mesh.cpp \
	mesh_io.cpp \
	MeshOpt.cpp \
	off_io.cpp \
	OneDefectRemeshing.cpp \
	QuadCleanUp.cpp \
	QuadPatches.cpp \
        QuadMesh.cpp \
	Singlet.cpp \
	SwapTriEdges.cpp \
	SwapQuadEdges.cpp \
	tfiblend.cpp \
	tri_io.cpp \
	circumcenter.cpp \
	GeomPredicates.cpp \
	MeshRefine2D.cpp \
	MeshQuality.cpp  \
	vtk_io.cpp      \
	xml_io.cpp     \
        Tri2Quad.cpp \
        cubit_io.cpp 
 
  

# The list of header files which are to be installed
nobase_libJaalQuadMesher_la_include_HEADERS = \
        meshkit/QuadMesh.hpp \
        meshkit/DijkstraShortestPath.hpp \
        meshkit/StopWatch.hpp   \
        meshkit/circumcenter.hpp \
        meshkit/GeomPredicates.hpp \
        meshkit/Tri2Quad.hpp \
        meshkit/tfiblend.hpp \
        meshkit/MeshRefine2D.hpp \
        meshkit/SwapEdges.hpp \
        meshkit/QuadCleanUp.hpp \
        meshkit/Mesh.hpp \
        meshkit/JaalMoabConverter.hpp \
        meshkit/DualGraph.hpp \
        meshkit/FastArea.hpp \
        meshkit/BinaryTree.hpp \
        meshkit/ObjectPool.hpp \
        meshkit/QuadPatchRemesh.hpp \
        meshkit/cycle.hpp \
        meshkit/basic_math.hpp \
        meshkit/iMeshRefine2D.hpp \
        meshkit/myany.hpp

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
