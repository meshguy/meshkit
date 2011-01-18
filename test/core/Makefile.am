SUBDIRS = . 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

LDADD = $(top_builddir)/src/libMeshKit.la

AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/core 

TESTS = test_graph

check_PROGRAMS = $(TESTS)

MOSTLYCLEANFILES = 

# Tests and such

test_graph_SOURCES = test_graph.cpp
test_graph_DEPENDENCIES = \
	$(LDADD)

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS)
AM_LDFLAGS += $(IGEOM_LTFLAGS)
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

