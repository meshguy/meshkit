# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libCoreGen.la
         
SUBDIRS = 

# The directory in which to install headers
libCoreGen_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libCoreGen_la_SOURCES = CoreGen.cpp \
        copymove_parallel.cpp
# The list of header files which are to be installed
nobase_libCoreGen_la_include_HEADERS = \
        meshkit/CoreGen.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../../core \
              -I$(srcdir)/../../lemon \
              -I$(srcdir)/../../utils \
              -I$(srcdir)/../../algs \
              -I$(top_builddir) \
              -I$(top_srcdir)/utils 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

LDADD = $(top_builddir)/src/libMeshKit.la
