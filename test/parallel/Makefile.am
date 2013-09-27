SUBDIRS = .

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/algs/CoreGen \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils \
               -I$(top_srcdir)/src/extern/CAMAL

BROKEN_TESTS = test_coregen

TESTS = test_coregen

if HAVE_CAMAL
TESTS += test_parallelmesher
endif

check_PROGRAMS = $(TESTS)

test_parallelmesher_SOURCES = test_parallelmesher.cpp
test_coregen_SOURCES = test_coregen.cpp


# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
