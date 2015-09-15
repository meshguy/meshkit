SUBDIRS = . 
MOSTLYCLEANFILES = 5bout.h5m AfterOneSetup.eps AfterTFISetup.eps AfterWinslow.vtk BeforeHex.vtk \
BeforeVolumeImprove.h5m InitialMapping.vtk OneToOneSwept.h5m SCDmesh.h5m SCDmesh1.vtk SCDmesh2.vtk \
SCDmesh3.vtk TFIMapping.vtk assygen_default.jou assygen_default.brep assygen_default.sat assygen_default.template.jou graph.eps \
intassign.h5m mesquite.gpt test.h5m volumesIce.h5m OneToOneSwept.vtk test_postbl.inp.log 5bout.exo ebmesh_out.vtk input.vtk \
assygen_default.stp cyl_mod.h5m cgd.h5m coregen_d.h5m CoreGen_default.log CoreGen_default.makefile

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data
if HAVE_TRIANGLE
AM_LDFLAGS += $(TRIANGLE_LIBS)
endif

if HAVE_TETGEN
AM_LDFLAGS +=$(TETGEN_LIBS)
endif

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/extern/CAMAL \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/algs/Sweep \
               -I$(top_srcdir)/src/algs/Qslim \
               -I$(top_srcdir)/src/algs/PostBL \
               -I$(top_srcdir)/src/algs/CoreGen \
               -I$(top_srcdir)/src/algs/AssyGen \
               -I$(top_srcdir)/src/algs/IntervalAssignment \
               -I$(top_srcdir)/src/extern/triangle \
                -I$(top_srcdir)/src/extern/tetgen \
     	       -I$(top_srcdir)/src/algs/make_watertight \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif


postbl_SOURCES = postbl.cpp

coregen_SOURCES = coregen.cpp

assygen_SOURCES = assygen.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

bin_PROGRAMS =


bin_PROGRAMS += assygen coregen postbl


