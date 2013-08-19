# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libPostBL.la
         
SUBDIRS = 

# The directory in which to install headers
libPostBL_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libPostBL_la_SOURCES = \
	mstream.hpp    PostBL.cpp    mstream.cpp  algorithm2.cpp
# The list of header files which are to be installed
nobase_libPostBL_la_include_HEADERS = \
        ../meshkit/PostBL.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../../core \
              -I$(srcdir)/../../lemon \
              -I$(srcdir)/../../utils \
              -I$(top_builddir) \
              -I$(top_srcdir)/utils 

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
