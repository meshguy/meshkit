SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la
DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/extern/mesquite \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir) \
               -I$(top_srcdir)/utils \
               $(MESQUITE_DEFINES) \
               $(MESQUITE_CPPFLAGS)
if HAVE_OCC_GEOMETRY_ENGINE
TESTS = test_mesquiteopt
check_PROGRAMS = $(TESTS)
endif
test_mesquiteopt_SOURCES = test_mesquiteopt.cpp
test_mesquiteopt_LDADD = ../../../src/libMeshKit.la

MOSTLYCLEANFILES = 

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check


