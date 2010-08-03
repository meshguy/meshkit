# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Tests and such
check_PROGRAMS = CopyMesh_test ExtrudeMesh_test EBMesh_test \
		 ProjectShell_test 
TESTS = $(check_PROGRAMS)

CopyMesh_test_SOURCES = CopyMesh.cpp
CopyMesh_test_DEPENDENCIES = ../CopyMesh.hpp ../libMKalgs.la
CopyMesh_test_CXXFLAGS = $(CXXFLAGS)
CopyMesh_test_LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp
ExtrudeMesh_test_DEPENDENCIES = ../ExtrudeMesh.hpp ../libMKalgs.la
ExtrudeMesh_test_CXXFLAGS = $(CXXFLAGS)
ExtrudeMesh_test_LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

EBMesh_test_SOURCES = EBMesh.cpp
EBMesh_test_DEPENDENCIES = ../EBMesh.hpp ../libMKalgs.la
EBMesh_test_CXXFLAGS = $(CXXFLAGS)
EBMesh_test_LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

ProjectShell_test_SOURCES = ProjectShell.cpp
ProjectShell_test_DEPENDENCIES = ../ProjectShell.hpp ../libMKalgs.la
ProjectShell_test_CXXFLAGS = $(CXXFLAGS)
ProjectShell_test_LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils -I$(top_srcdir)/algs

if HAVE_IGEOM
INCLUDES += $(IGEOM_INCLUDES)
CPPFLAGS += $(IGEOM_CPPFLAGS) -DCGM
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
INCLUDES += $(IREL_INCLUDES)
CPPFLAGS += $(IREL_CPPFLAGS)
endif

if HAVE_MOAB
CPPFLAGS += -DMOAB
endif