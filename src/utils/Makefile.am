# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitutils.la
SUBDIRS = 

# The directory in which to install headers
libMeshKitutils_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitutils_la_SOURCES = \
	iMesh.cpp

# The list of header files which are to be installed
nobase_libMeshKitutils_la_include_HEADERS = \
	meshkit/iGeom.hh                 \
	meshkit/iMesh.hh                 \
	meshkit/iRel.hh                  \
	meshkit/RegisterMeshOp.hpp         


# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core -I$(srcdir)/../utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/ -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
