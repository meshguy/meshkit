# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = Qslim test

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKalgs.la

# The directory in which to install headers
libMKalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKalgs_la_SOURCES = \
	CopyMesh.cpp \
	MoveMesh.cpp \
	ExtrudeMesh.cpp \
	MergeMesh.cpp \
	Transform.cpp \
	vec_utils.cpp \
	vec_utils.hpp \
	EBMesh.cpp \
	ProjectShell.cpp \
	EdgeMesher.cpp \
	OneToOneSwept.cpp \
	CESets.cpp \
	CopyUtils.cpp \
	CopyUtils.hpp

# The list of header files which are to be installed
libMKalgs_la_include_HEADERS = \
	CopyMesh.hpp \
	MoveMesh.hpp \
	ExtrudeMesh.hpp \
	MergeMesh.hpp \
	Transform.hpp \
	EBMesh.hpp \
	ProjectShell.hpp \
	EdgeMesher.hpp \
	OneToOneSwept.hpp \
	CESets.hpp

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils

if HAVE_IGEOM
INCLUDES += $(IGEOM_INCLUDES)
CPPFLAGS += $(IGEOM_CPPFLAGS) -DCGM
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
INCLUDES += $(IREL_INCLUDES)
CPPFLAGS += $(IREL_CPPFLAGS)
endif

if HAVE_MOAB
CPPFLAGS += -DMOAB
endif

