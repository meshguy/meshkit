# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS =  QuadMesher Qslim AssyGen CoreGen Sweep PostBL

# Things to build
noinst_LTLIBRARIES = libMeshKitalgs.la
libMeshKitalgs_la_LIBADD = Qslim/libQslim.la \
                           QuadMesher/libJaalQuadMesher.la \
                           Sweep/libSweep.la \
                           AssyGen/libAssyGen.la \
			   CoreGen/libCoreGen.la \
                           PostBL/libPostBL.la                        

if ENABLE_intassign
  SUBDIRS += IntervalAssignment
  AM_CPPFLAGS += -I$(top_srcdir)/src/algs/IntervalAssignment $(IPOPT_CPPFLAGS)
  libMeshKitalgs_la_LIBADD += IntervalAssignment/libIntervalAssignment.la
  AM_LDFLAGS += $(IPOPT_LDFLAGS)
  LIBS += $(IPOPT_LIBS)
endif

# The directory in which to install headers
libMeshKitalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitalgs_la_SOURCES = \
	register_algs.cpp \
	EdgeMesher.cpp \
	VertexMesher.cpp \
	SCDMesh.cpp \
	EBMesher.cpp \
	CopyMesh.cpp \
	CopyGeom.cpp \
	MeshOpTemplate.cpp \
	MergeMesh.cpp \
	ExtrudeMesh.cpp \
	CESets.cpp \
        Transform.cpp
if HAVE_FBIGEOM
libMeshKitalgs_la_SOURCES += \
  MBGeomOp.cpp \
  MBSplitOp.cpp \
  MBVolOp.cpp
endif

# The list of header files which are to be installed
nobase_libMeshKitalgs_la_include_HEADERS = \
	meshkit/EdgeMesher.hpp \
	meshkit/VertexMesher.hpp \
	meshkit/SCDMesh.hpp \
	meshkit/EBMesher.hpp \
	meshkit/CopyMesh.hpp \
	meshkit/CopyGeom.hpp \
	meshkit/MeshOpTemplate.hpp \
	meshkit/MergeMesh.hpp \
	meshkit/ExtrudeMesh.hpp \
	meshkit/CESets.hpp \
	meshkit/TransformBase.hpp \
        meshkit/Transform.hpp
if HAVE_FBIGEOM
nobase_libMeshKitalgs_la_include_HEADERS += \
	meshkit/MBGeomOp.hpp \
	meshkit/MBSplitOp.hpp \
	meshkit/MBVolOp.hpp
endif

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) \
              -I$(srcdir)/../core \
              -I$(srcdir)/../utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/ \
              -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon \
              -I$(srcdir)/QuadMesher \
              -I$(srcdir)/Sweep \
              -I$(srcdir)/AssyGen \
	      -I$(srcdir)/CoreGen \
	      -I$(srcdir)/PostBL \
	      -I$(srcdir)/Qslim

if PARALLEL
if HAVE_PARALLEL_MOAB
if HAVE_PARALLEL_CGM
  libMeshKitalgs_la_SOURCES += ParallelMesher.cpp ParExchangeMesh.cpp ParSendPostSurfMesh.cpp ParRecvSurfMesh.cpp
  nobase_libMeshKitalgs_la_include_HEADERS += meshkit/ParallelMesher.hpp meshkit/ParExchangeMesh.hpp meshkit/ParSendPostSurfMesh.hpp meshkit/ParRecvSurfMesh.hpp
  AM_CPPFLAGS += -I$(srcdir)/../../extern/CAMAL
endif
endif
endif
