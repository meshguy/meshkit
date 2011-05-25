SUBDIRS = core algs extern

if PARALLEL
if HAVE_PARALLEL_MOAB
if HAVE_PARALLEL_CGM
SUBDIRS += parallel
endif
endif
endif

build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

TESTS = test_testutil
check_PROGRAMS = $(TESTS)

test_testutil_SOURCES = TestUtil.hpp test_testutil.cpp
