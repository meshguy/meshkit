SUBDIRS = . 

LDADD = $(top_builddir)/src/libMeshKit.la

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/test \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/utils \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/algs \
               -I$(top_srcdir)/src/lemon -I$(top_srcdir)

TESTS = test_edgemesh test_ebmesh

check_PROGRAMS = $(TESTS) test_copymesh test_scdmesh

MOSTLYCLEANFILES = 

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

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

