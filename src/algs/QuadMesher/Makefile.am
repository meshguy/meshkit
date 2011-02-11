# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitalgs.la
         
SUBDIRS = 

# The directory in which to install headers
libMeshKitalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitalgs_la_SOURCES = \
	AdfrontClean.cc  \
	BinaryTree.cc \
	Bridges.cc \
	DijkstraShortestPath.cc \
	Dimonds.cc \
	Doublet.cc \
	DualGraph.cc \
	edmond_matching.cc \
	FaceClose.cc \
	laplace.cc \
	matching.cc \
	Mesh2Moab.cc \
	Mesh.cc \
	meshclean.cc \
	mesh_io.cc \
	MeshOpt.cc \
	meshViewer.cc \
	mreport.cc \
	off_io.cc \
	OneDefectRemeshing.cc \
	patterns.cc \
	QuadCleanUp.cc \
	QuadPatches.cc \
	Singlet.cc \
	SwapQuadEdge.cc \
	SwapTriEdge.cc \
	tfiblend.cc \
	Tri2Quad.cc \
	tri_io.cc \
	Tunnels.cc \
	vtk_io.cc \
	circumcenter.cpp \
	geomPredicates.cpp \
	MeshRefine2D.cpp \
	meshuntangle.cpp \
	predicates.cpp \
	QuadQuality.cpp \

# The list of header files which are to be installed
nobase_libMeshKitalgs_la_include_HEADERS = \
	meshkit/basic_math.h
	meshkit/BinaryTree.h
	meshkit/circumcenter.h
	meshkit/cycle.h
	meshkit/DijkstraShortestPath.h
	meshkit/DualGraph.h
	meshkit/FastArea.h
	meshkit/geomPredicates.h
	meshkit/iMeshRefine2D.h
	meshkit/Mesh.h
	meshkit/MeshRefine2D.h
	meshkit/ObjectPool.h
	meshkit/QuadCleanUp.h
	meshkit/SwapTriEdge.h
	meshkit/tfiblend.h
	meshkit/Tri2Quad.h

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/ -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
