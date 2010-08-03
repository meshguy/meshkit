# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Tests and such
check_PROGRAMS = \
	CopyMesh_test \
	EBMesh_test \
	ExtrudeMesh_test \
	ProjectShell_test 
TESTS = $(check_PROGRAMS)

LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

CopyMesh_test_SOURCES = CopyMesh.cpp TestFramework.hpp
EBMesh_test_SOURCES = EBMesh.cpp TestFramework.hpp
ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp TestFramework.hpp
ProjectShell_test_SOURCES = ProjectShell.cpp TestFramework.hpp

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