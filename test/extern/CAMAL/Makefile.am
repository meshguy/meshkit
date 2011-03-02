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

MOSTLYCLEANFILES = 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# Tests and such
if HAVE_CAMAL_TRIADV
test_triadvance_SOURCES = test_triadvance.cpp
test_triadvance_DEPENDENCIES = \
	$(LDADD)
test_mbgeom_SOURCES = test_mbgeom.cpp
test_mbgeom_DEPENDENCIES = $(LDADD)
TESTS += test_triadvance test_mbgeom
endif

if HAVE_CAMAL_PAVER
test_paver_SOURCES = test_paver.cpp
test_paver_DEPENDENCIES = \
	$(LDADD)
TESTS += test_paver
endif

if HAVE_CAMAL_TET
test_tetmesher_SOURCES = test_tetmesher.cpp
test_tetmesher_DEPENDENCIES = \
	$(LDADD)
TESTS += test_tetmesher
endif

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


