# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = \
	3k-tri-cube.h5m \
	assygen_default.cub \
	m1.vtk \
	m2.vtk \
	sphere.stp 

MOSTLYCLEANFILES = \
	ProjShell_3k-tri-cube.h5m \
	t1.h5m \
	t2.h5m \
	t3.h5m \
	t4.h5m \
	t5.h5m \
	intx.vtk

# Tests and such
check_PROGRAMS = \
	CopyMesh_test \
	ExtrudeMesh_test \
	ProjectShell_test \
    	IntersectMesh_test \
    	CopyGeom_test AddZTag
TESTS = $(check_PROGRAMS)

LDADD = ../libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

CopyMesh_test_SOURCES = CopyMesh.cpp TestFramework.hpp 
ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp TestFramework.hpp
ProjectShell_test_SOURCES = ProjectShell.cpp TestFramework.hpp
IntersectMesh_test_SOURCES = IntersectMesh_test.cpp TestFramework.hpp
CopyGeom_test_SOURCES = CopyGeom_test.cpp TestFramework.hpp
pointset_SOURCES = pointset.cpp
check_PROGRAMS += pointset
AddZTag_SOURCES = kp.AddZTag_test.cpp

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
