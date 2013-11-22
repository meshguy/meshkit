# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
lib_LTLIBRARIES = libMKplibs.la

# The directory in which to install headers
libMKplibs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKplibs_la_SOURCES = \
	EdgeMesher.cc \
	SurfMesher.cc \
	VolMesher.cc

# The list of header files which are to be installed
libMKplibs_la_include_HEADERS = \
	Mesh.h \
	itaps.h \
	EdgeMesher.h \
	SurfMesher.h \
	VolMesher.h

# Tests and such
check_PROGRAMS = TetGen_test
TESTS = $(check_PROGRAMS)

TetGen_test_SOURCES = TetGen_test.cpp
TetGen_test_DEPENDENCIES = libMKplibs.la
TetGen_test_CXXFLAGS = -DTEST $(CXXFLAGS)
TetGen_test_LDADD = libMKplibs.la $(TETGEN_LIBS) $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) -I$(top_srcdir) -I$(top_builddir)

if HAVE_TETGEN
CPPFLAGS += $(TETGEN_CPPFLAGS) -DTETLIBRARY
#LDFLAGS += $(TETGEN_LDFLAGS)
endif
