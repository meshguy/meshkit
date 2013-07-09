SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/example \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/extern/CAMAL \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/algs/Sweep \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils

BROKEN_EXAMPLES = 

EXAMPLES = example_basic

if ENABLE_intassign
EXAMPLES += test_intassign
endif

check_PROGRAMS = $(EXAMPLES) $(BROKEN_EXAMPLES)

MOSTLYCLEANFILES = 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# Tests and such
example_basic_SOURCES = example_basic.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'EXAMPLES_ENVIRONMENT=: ' check


