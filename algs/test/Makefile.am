# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Tests and such
check_PROGRAMS = \
	CopyMesh_test \
	ExtrudeMesh_test \
	ProjectShell_test \
    	IntersectMesh_test \
    	CopyGeom_test
TESTS = $(check_PROGRAMS)

LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

CopyMesh_test_SOURCES = CopyMesh.cpp TestFramework.hpp 
#EBMesh_test_SOURCES = EBMesh.cpp TestFramework.hpp
ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp TestFramework.hpp
ProjectShell_test_SOURCES = ProjectShell.cpp TestFramework.hpp
IntersectMesh_test_SOURCES = IntersectMesh_test.cpp TestFramework.hpp
CopyGeom_test_SOURCES = CopyGeom_test.cpp TestFramework.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir)\
               -I$(top_srcdir) \
               -I$(top_builddir) \
               -I$(top_srcdir)/utils \
               -I$(top_srcdir)/algs

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

if HAVE_TRIANGLE
BedJako_SOURCES = BedJako.cpp TestFramework.hpp
check_PROGRAMS += BedJako
AM_CPPFLAGS += $(TRIANGLE_CPPFLAGS)
LDADD += $(TRIANGLE_LIBS)
endif
