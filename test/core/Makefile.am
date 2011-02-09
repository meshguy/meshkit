SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir)

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

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

