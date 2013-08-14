SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/utils \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir)

TESTS = test_graph list_meshops

check_PROGRAMS = $(TESTS)

MOSTLYCLEANFILES = graph.eps

# Tests and such

test_graph_SOURCES = test_graph.cpp

list_meshops_SOURCES = list_meshops.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

