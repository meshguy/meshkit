# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libAssyMesher.la
         
SUBDIRS = 

# The directory in which to install headers
libAssyMesher_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libAssyMesher_la_SOURCES = AssyMesher.cpp
# The list of header files which are to be installed
nobase_libAssyMesher_la_include_HEADERS = \
        meshkit/AssyMesher.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core \
              -I$(top_srcdir)/src/utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/src/lemon \
              -I$(top_srcdir)/src/algs/Sweep

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

LDADD = $(top_builddir)/src/libMeshKit.la
