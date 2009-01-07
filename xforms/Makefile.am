# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKxforms.la

# The directory in which to install headers
libMKxforms_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKxforms_la_SOURCES = 

# The list of header files which are to be installed
libMKxforms_la_include_HEADERS = 

# Tests and such
check_PROGRAMS = mk_xforms_test
TESTS = $(check_PROGRAMS)

mk_xforms_test_SOURCES = mk_xforms_test.cpp
mk_xforms_test_LDADD = libMKxforms.la $(iMesh_LIBS) $(iGeom_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
INCLUDES += $(iGeom_INCLUDES)
CPPFLAGS += $(iGeom_CPPFLAGS)
CXXFLAGS += $(iGeom_CXXFLAGS)
CFLAGS += $(iGeom_CFLAGS)
LDFLAGS += $(iGeom_LDFLAGS)
endif

if HAVE_IMESH
INCLUDES += $(iMesh_INCLUDES)
CPPFLAGS += $(iMesh_CPPFLAGS)
CXXFLAGS += $(iMesh_CXXFLAGS)
CFLAGS += $(iMesh_CFLAGS)
LDFLAGS += $(iMesh_LDFLAGS)
endif

