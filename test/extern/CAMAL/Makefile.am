SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la $(top_builddir)/src/extern/libMeshKitextern.la

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

@IGEOM_CONFIG_OPTIONS_INCLUDE@
@IMESH_CONFIG_OPTIONS_INCLUDE@
@IREL_CONFIG_OPTIONS_INCLUDE@

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS)
LIBS += $(IGEOM_LIBS)
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES) $(IMESH_CPPFLAGS)
LIBS += $(IMESH_LIBS)
endif

if HAVE_IREL
AM_CPPFLAGS += $(IREL_INCLUDES) $(IREL_CPPFLAGS)
LIBS += $(IREL_LIBS)
endif

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

if HAVE_CAMAL_TET
test_tetmesher_SOURCES = test_tetmesher.cpp
test_tetmesher_DEPENDENCIES = \
	$(LDADD)
TESTS += test_tetmesher
endif

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


