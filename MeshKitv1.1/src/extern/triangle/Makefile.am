# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKittriangle.la
SUBDIRS = 

# The directory in which to install headers
libMeshKittriangle_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKittriangle_la_SOURCES = TriangleMesher.cpp

# The list of header files which are to be installed
nobase_libMeshKittriangle_la_include_HEADERS = meshkit/TriangleMesher.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core -I$(top_srcdir)/src/utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/utils -I$(top_srcdir)/src/extern/triangle \
              -I$(top_srcdir)/src/lemon

if HAVE_TRIANGLE
AM_CPPFLAGS += $(TRIANGLE_CPPFLAGS)
LDFLAGS += $(TRIANGLE_LDFLAGS)
endif

