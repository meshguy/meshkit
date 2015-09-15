# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKittetgen.la
SUBDIRS = 

# The directory in which to install headers
libMeshKittetgen_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKittetgen_la_SOURCES = TGTetMesher.cpp

# The list of header files which are to be installed
nobase_libMeshKittetgen_la_include_HEADERS = meshkit/TGTetMesher.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core -I$(top_srcdir)/src/utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/utils -I$(top_srcdir)/src/extern/tetgen \
              -I$(top_srcdir)/src/lemon

if HAVE_TETGEN
AM_CPPFLAGS += $(TETGEN_CPPFLAGS)
LDFLAGS += $(TETGEN_LDFLAGS)
endif

