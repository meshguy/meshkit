SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/utils \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) 

TESTS = test_graph list_meshops test_load_grps test_get_entities_by_dim test_imesh

check_PROGRAMS = $(TESTS)

MOSTLYCLEANFILES = graph.eps

# Tests and such

test_imesh_SOURCES = test_imesh.cpp

test_graph_SOURCES = test_graph.cpp

list_meshops_SOURCES = list_meshops.cpp

test_load_grps_SOURCES = test_load_grps.cpp

test_get_entities_by_dim_SOURCES = test_get_entities_by_dim.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

