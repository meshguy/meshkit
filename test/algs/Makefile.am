SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

TESTS = test_edgemesh test_ebmesh test_scdmesh

check_PROGRAMS = $(TESTS) test_copymesh

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

test_edgemesh_SOURCES = test_edgemesh.cpp
test_edgemesh_DEPENDENCIES = \
	$(LDADD)

test_ebmesh_SOURCES = test_ebmesh.cpp
test_ebmesh_DEPENDENCIES = \
	$(LDADD)

test_copymesh_SOURCES = test_copymesh.cpp
test_copymesh_DEPENDENCIES = \
	$(LDADD)

test_scdmesh_SOURCES = test_scdmesh.cpp
test_scdmesh_DEPENDENCIES = \
	$(LDADD)

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


