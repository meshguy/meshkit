SUBDIRS = .
MOSTLYCLEANFILES =

LDADD = $(top_builddir)/src/libMeshKit.la

# Define data directory and pass its definition to the preprocessor
DATADIR = $(top_srcdir)/data
AM_CPPFLAGS += -DDATADIR=$(DATADIR)

AM_CPPFLAGS += -I$(top_srcdir)/src/algs/AdvFront \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/lemon \
               -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/test


TESTS = freeZone ruleAppObj2D triRule planeProj 

check_PROGRAMS = $(TESTS)

freeZone_SOURCES = freeZone.cpp
planeProj_SOURCES = planeProj.cpp
ruleAppObj2D_SOURCES = ruleAppObj2D.cpp
triRule_SOURCES = triRule.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
