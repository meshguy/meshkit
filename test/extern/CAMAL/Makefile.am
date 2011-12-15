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
TESTS += test_triadvance
endif

if HAVE_CAMAL_PAVER
test_paver_SOURCES = test_paver.cpp
test_paver_DEPENDENCIES = \
	$(LDADD)
TESTS += test_paver
if HAVE_FBIGEOM
test_mbgeom_SOURCES = test_mbgeom.cpp
test_mbgeom_DEPENDENCIES = $(LDADD)
test_mbgeomvar_SOURCES = test_mbgeomvar.cpp
test_mbgeomvar_DEPENDENCIES = $(LDADD)
test_mb2_SOURCES= test_mb2.cpp
test_mb2_DEPENDENCIES = $(LDADD)
TESTS += test_mbgeom test_mbgeomvar test_mb2
# test_mb2 is broken for the time being
endif
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


