# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS =

# Things to build
noinst_LTLIBRARIES = libMeshKitmesquite.la
libMeshKitmesquite_la_LIBADD = $(MESQUITE_LIBS)

# The directory in which to install headers
libMeshKitmesquite_la_includedir = $(includedir)

libMeshKitmesquite_la_SOURCES = \
	MesquiteOpt.cpp \
	FreeSmoothDomain.cpp \
	FreeSmoothDomain.hpp

# The list of header files which are to be installed
nobase_libMeshKitmesquite_la_include_HEADERS = \
	meshkit/MesquiteOpt.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core -I$(top_srcdir)/src/utils \
              -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils \
              -I$(top_srcdir)/src/extern/mesquite \
              -I$(top_srcdir)/src/lemon $(MESQUITE_DEFINES) \
              $(MESQUITE_CPPFLAGS) 

AM_LDFLAGS += $(MESQUITE_LDFLAGS)

