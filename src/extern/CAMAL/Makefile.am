# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitCAMAL.la
SUBDIRS = 

# The directory in which to install headers
libMeshKitCAMAL_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitCAMAL_la_SOURCES = \
	CAMALGeomEval.cpp \
	CAMALGeomEval.hpp

# The list of header files which are to be installed
nobase_libMeshKitCAMAL_la_include_HEADERS = 

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core -I$(top_srcdir)/src/utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/utils -I$(top_srcdir)/src/extern/CAMAL \
              -I$(top_srcdir)/src/lemon

if HAVE_CAMAL
AM_CPPFLAGS += -DCAMAL_VERSION=$(CAMAL_VERSION) $(CAMAL_CPPFLAGS)
LDFLAGS += $(CAMAL_LDFLAGS)
endif

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
