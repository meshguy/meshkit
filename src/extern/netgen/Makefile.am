# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitnetgen.la
SUBDIRS =

# The directory in which to install headers
libMeshKitnetgen_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitnetgen_la_SOURCES = NGTetMesher.cpp

# The list of header files which are to be installed
nobase_libMeshKitnetgen_la_include_HEADERS = meshkit/NGTetMesher.hpp

if HAVE_OCC_NETGEN
libMeshKitnetgen_la_SOURCES += NGTriMesher.cpp
nobase_libMeshKitnetgen_la_include_HEADERS += meshkit/NGTriMesher.hpp
endif

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core -I$(top_srcdir)/src/utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/utils -I$(top_srcdir)/src/extern/netgen \
              -I$(top_srcdir)/src/lemon

if HAVE_NETGEN
AM_CPPFLAGS += $(NETGEN_CPPFLAGS)
LDFLAGS += $(NETGEN_LDFLAGS)
endif
