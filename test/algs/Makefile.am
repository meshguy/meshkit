SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/extern/CAMAL \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

BROKEN_TESTS = test_jaalquad test_assygen

TESTS = test_edgemesh \
        test_tfimapping \
        test_onetooneswept \
        test_scdmesh \
        test_copymesh \
        test_copygeom \
        test_extrudemesh \
        test_qslim \
        test_fbgeom \
        test_ebmesh \
        test_assygen \
        test_TriangleMesher \
        test_MBGeomOp \
        test_MBSplitOp

check_PROGRAMS = $(TESTS) $(BROKEN_TESTS)

MOSTLYCLEANFILES = SCDmesh1.vtk SCDmesh2.vtk SCDmesh3.vtk SCDmesh.h5m

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# Tests and such

test_edgemesh_SOURCES = test_edgemesh.cpp
test_tfimapping_SOURCES = test_tfimapping.cpp
test_onetooneswept_SOURCES = test_onetooneswept.cpp
test_ebmesh_SOURCES = test_ebmesh.cpp
test_copymesh_SOURCES = test_copymesh.cpp
test_copygeom_SOURCES = test_copygeom.cpp
test_extrudemesh_SOURCES = test_extrudemesh.cpp
test_scdmesh_SOURCES = test_scdmesh.cpp
test_qslim_SOURCES = test_qslim.cpp
test_jaalquad_SOURCES = test_jaalquad.cpp
test_assygen_SOURCES = test_assygen.cpp
test_TriangleMesher_SOURCES = test_TriangleMesher.cpp
test_MBGeomOp_SOURCES = test_MBGeomOp.cpp
test_MBSplitOp_SOURCES = test_MBSplitOp.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


