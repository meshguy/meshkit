# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

SUBDIRS = Qslim test make_watertight

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKalgs.la

#libMKalgs_la_LIBADD = make_watertight/libMakeWatertight.la

# The directory in which to install headers
libMKalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKalgs_la_SOURCES = \
	CopyMesh.cpp \
	MoveMesh.cpp \
	ExtrudeMesh.cpp \
	MergeMesh.cpp \
	Transform.cpp \
	vec_utils.cpp \
	vec_utils.hpp \
	EBMesh.cpp \
	ProjectShell.cpp \
	IntersectMesh.cpp \
	EdgeMesher.cpp \
	OneToOneSwept.cpp \
	CESets.cpp \
	CopyUtils.cpp \
	CopyUtils.hpp

# The list of header files which are to be installed
libMKalgs_la_include_HEADERS = \
	CopyMesh.hpp \
	MoveMesh.hpp \
	ExtrudeMesh.hpp \
	MergeMesh.hpp \
	Transform.hpp \
	EBMesh.hpp \
	ProjectShell.hpp \
	IntersectMesh.hpp \
	EdgeMesher.hpp \
	OneToOneSwept.hpp \
	CESets.hpp

# Tests and such
check_PROGRAMS = CopyMesh_test ExtrudeMesh_test EBMesh_test \
		 ProjectShell_test matching
TESTS = $(check_PROGRAMS)

matching_SOURCES = QuadMesh/matching.cc
matching_DEPENDENCIES = libMKalgs.la
matching_CXXFLAGS = $(CXXFLAGS)
matching_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS)

CopyMesh_test_SOURCES = CopyMesh.cpp
CopyMesh_test_DEPENDENCIES = CopyMesh.hpp libMKalgs.la
CopyMesh_test_CXXFLAGS = -DTEST $(CXXFLAGS)
CopyMesh_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp
ExtrudeMesh_test_DEPENDENCIES = ExtrudeMesh.hpp libMKalgs.la
ExtrudeMesh_test_CXXFLAGS = -DTEST $(CXXFLAGS)
ExtrudeMesh_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

EBMesh_test_SOURCES = EBMesh_test.cpp
EBMesh_test_DEPENDENCIES = EBMesh.hpp libMKalgs.la
EBMesh_test_CXXFLAGS = -DTEST $(CXXFLAGS)
EBMesh_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

ProjectShell_test_SOURCES = ProjectShell_test.cpp
ProjectShell_test_DEPENDENCIES = ProjectShell.hpp libMKalgs.la
ProjectShell_test_CXXFLAGS = -DTEST $(CXXFLAGS)
ProjectShell_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS) -DCGM
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES) $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
AM_CPPFLAGS += $(IREL_INCLUDES) $(IREL_CPPFLAGS)
endif

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

