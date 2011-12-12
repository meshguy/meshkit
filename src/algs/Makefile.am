# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS =  QuadMesher Qslim AssyGen Triangle Sweep

# Things to build
noinst_LTLIBRARIES = libMeshKitalgs.la
libMeshKitalgs_la_LIBADD = Qslim/libMKQslim.la \
                           QuadMesher/libJaalQuadMesher.la \
                           Sweep/libSweep.la \
                           AssyGen/libAssyGen.la \
                           Triangle/libMKTriangle.la

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
    MergeMesh.cpp \
	ExtrudeMesh.cpp \
	CESets.cpp \
	Transform.cpp \
	CopyUtils.cpp \
	CopyUtils.hpp 
if HAVE_FBIGEOM
libMeshKitalgs_la_SOURCES += \
  MBGeomOp.cpp \
  MBSplitOp.cpp
endif

# The list of header files which are to be installed
nobase_libMeshKitalgs_la_include_HEADERS = \
	meshkit/EdgeMesher.hpp \
	meshkit/VertexMesher.hpp \
	meshkit/SCDMesh.hpp \
	meshkit/EBMesher.hpp \
	meshkit/CopyMesh.hpp \
	meshkit/CopyGeom.hpp \
	meshkit/MergeMesh.hpp \
	meshkit/ExtrudeMesh.hpp \
	meshkit/CESets.hpp \
	meshkit/TransformBase.hpp \
	meshkit/Transform.hpp 
if HAVE_FBIGEOM
nobase_libMeshKitalgs_la_include_HEADERS += \
	meshkit/MBGeomOp.hpp \
	meshkit/MBSplitOp.hpp
endif

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core \
              -I$(srcdir)/../utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/ \
              -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon \
              -I$(srcdir)/QuadMesher \
              -I$(srcdir)/Sweep \
	      -I$(srcdir)/AssyGen


if PARALLEL
if HAVE_PARALLEL_MOAB
if HAVE_PARALLEL_CGM
  libMeshKitalgs_la_SOURCES += ParallelMesher.cpp ParExchangeMesh.cpp
  nobase_libMeshKitalgs_la_include_HEADERS += meshkit/ParallelMesher.hpp meshkit/ParExchangeMesh.hpp
  AM_CPPFLAGS += -I$(srcdir)/../../extern/CAMAL
endif
endif
endif
