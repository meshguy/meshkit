# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libPbl.la
         
SUBDIRS = 

# The directory in which to install headers
libPbl_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libPbl_la_SOURCES =  \
        ./Pbl.cpp
# The list of header files which are to be installed
nobase_libPbl_la_include_HEADERS = \
        ../meshkit/Pbl.hpp

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
