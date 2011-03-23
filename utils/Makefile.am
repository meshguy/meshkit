# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
lib_LTLIBRARIES = libMKutils.la

# The directory in which to install headers
libMKutils_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKutils_la_SOURCES = \
	LocalTag.hpp \
	LocalSet.hpp \
	MKUtils.cpp

# The list of header files which are to be installed
libMKutils_la_include_HEADERS = \
	MKException.hpp \
	ArrayManager.hpp \
	SimpleArray.hpp \
	MKUtils.hpp 

# Tests and such
check_PROGRAMS = MKUtils_test
TESTS = $(check_PROGRAMS)

MKUtils_test_SOURCES = MKUtils_test.cpp
MKUtils_test_DEPENDENCIES = MKUtils.hpp libMKutils.la
MKUtils_test_CXXFLAGS = -DTEST -DTEMPLATE_DEFS_INCLUDED $(CXXFLAGS)
MKUtils_test_LDADD = libMKutils.la $(IMESH_LIBS) $(IGEOM_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) -I$(top_srcdir) -I$(top_builddir)

