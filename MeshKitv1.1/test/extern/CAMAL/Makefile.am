SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la
DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/extern/CAMAL \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

TESTS = 

check_PROGRAMS = $(TESTS)

MOSTLYCLEANFILES = graph.eps simpletet.h5m threeholecube.h5m

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# Tests and such
if HAVE_CAMAL_TRIADV
test_camaltriadvance_SOURCES = test_camaltriadvance.cpp
test_camaltriadvance_DEPENDENCIES = \
	$(LDADD)
TESTS += test_camaltriadvance
endif

if HAVE_CAMAL_PAVER
test_camalpaver_SOURCES = test_camalpaver.cpp
test_camalpaver_DEPENDENCIES = \
	$(LDADD)
TESTS += test_camalpaver
if HAVE_FBIGEOM
test_camalmbgeom_SOURCES = test_camalmbgeom.cpp
test_camalmbgeom_DEPENDENCIES = $(LDADD)
test_camalmbgeomvar_SOURCES = test_camalmbgeomvar.cpp
test_camalmbgeomvar_DEPENDENCIES = $(LDADD)
test_camalmb2_SOURCES= test_camalmb2.cpp
test_camalmb2_DEPENDENCIES = $(LDADD)
TESTS += test_camalmbgeom test_camalmbgeomvar test_camalmb2
# test_camalmb2 is broken for the time being
endif
endif

if HAVE_CAMAL_TET
test_camaltetmesher_SOURCES = test_camaltetmesher.cpp
test_camaltetmesher_DEPENDENCIES = \
	$(LDADD)
TESTS += test_camaltetmesher
endif

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


