# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKutils.la

# The directory in which to install headers
libMKutils_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKutils_la_SOURCES = \
	MKUtils.cpp

# The list of header files which are to be installed
libMKutils_la_include_HEADERS = \
	ArrayManager.hpp \
	LocalTag.hpp \
	SimpleArray.hpp \
	MKUtils.hpp 

# Tests and such
check_PROGRAMS = MKUtils_test
TESTS = $(check_PROGRAMS)

MKUtils_test_SOURCES = MKUtils_test.cpp
MKUtils_test_DEPENDENCIES = MKUtils.hpp libMKutils.la
MKUtils_test_CXXFLAGS = -DTEST $(CXXFLAGS)
MKUtils_test_LDADD = libMKutils.la $(IMESH_LIBS) $(IGEOM_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS = $(DEFINES) -DSRCDIR=$(srcdir) -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
AM_CPPFLAGS += $(iGeom_INCLUDES) $(iGeom_CPPFLAGS)
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES) $(IMESH_CPPFLAGS)
endif

