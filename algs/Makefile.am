# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKalgs.la

# The directory in which to install headers
libMKalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKalgs_la_SOURCES = \
	CopyMesh.cpp \
	CopyVerts.cpp \
	ExtrudeMesh.cpp \
	MergeMesh.cpp \
	vec_utils.cpp \
	vec_utils.hpp

# The list of header files which are to be installed
libMKalgs_la_include_HEADERS = \
	CopyMesh.hpp \
	CopyVerts.hpp \
	ExtrudeMesh.hpp \
	MergeMesh.hpp

# Tests and such
check_PROGRAMS = CopyMesh_test ExtrudeMesh_test
TESTS = $(check_PROGRAMS)

CopyMesh_test_SOURCES = CopyMesh.cpp
CopyMesh_test_DEPENDENCIES = CopyMesh.hpp libMKalgs.la
CopyMesh_test_CXXFLAGS = -DTEST $(CXXFLAGS)
CopyMesh_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) 

ExtrudeMesh_test_SOURCES = ExtrudeMesh.cpp
ExtrudeMesh_test_DEPENDENCIES = ExtrudeMesh.hpp libMKalgs.la
ExtrudeMesh_test_CXXFLAGS = -DTEST $(CXXFLAGS)
ExtrudeMesh_test_LDADD = libMKalgs.la $(IMESH_LIBS) $(IGEOM_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
INCLUDES += $(iGeom_INCLUDES)
CPPFLAGS += $(iGeom_CPPFLAGS)
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
#CXXFLAGS += $(IMESH_CXXFLAGS)
#CFLAGS += $(IMESH_CFLAGS)
#LDFLAGS += $(IMESH_LDFLAGS)
endif

