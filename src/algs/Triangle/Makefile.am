# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMKTriangle.la

# The directory in which to install headers
libMKTriangle_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKTriangle_la_SOURCES = \
  TriangleMesher.cpp triangle.c triangle.h

# The list of header files which are to be installed (from algs/meshkit)
nobase_libMKTriangle_la_include_HEADERS = \
  ../meshkit/TriangleMesher.hpp 


# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../../core -I$(srcdir)/../../utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/ -I$(top_srcdir)/utils \
              -I$(srcdir)/../../lemon \
              -I$(srcdir)/..
   
# in order to build Triangle as a library           
AM_CPPFLAGS += -DTRILIBRARY
