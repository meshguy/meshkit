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


TESTS = pntTrnsfrmLnr2D freeZone ruleAppObj2D planeProj triRule dfltTriRules frontObj2D basicAlg2D triMeshOp

check_PROGRAMS = $(TESTS)

basicAlg2D_SOURCES = basicAlg2D.cpp
dfltTriRules_SOURCES = dfltTriRules.cpp
freeZone_SOURCES = freeZone.cpp
frontObj2D_SOURCES = frontObj2D.cpp
pntTrnsfrmLnr2D_SOURCES = pntTrnsfrmLnr2D.cpp
planeProj_SOURCES = planeProj.cpp
ruleAppObj2D_SOURCES = ruleAppObj2D.cpp
triMeshOp_SOURCES = triMeshOp.cpp
triRule_SOURCES = triRule.cpp

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
