SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la $(top_builddir)/src/extern/libMeshKitextern.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/extern/netgen \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

TESTS = 

check_PROGRAMS = $(TESTS)

MOSTLYCLEANFILES = 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# Tests and such
test_ngtetmesher_SOURCES = test_ngtetmesher.cpp
test_ngtetmesher_DEPENDENCIES = \
	$(LDADD)
TESTS += test_ngtetmesher

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


