# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

SUBDIRS = test make_watertight
#make_watertight

# Things to build
lib_LTLIBRARIES = libMKalgs.la

#libMKalgs_la_LIBADD = make_watertight/libMakeWatertight.la

# The directory in which to install headers
libMKalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKalgs_la_SOURCES = \
	CopyGeom.cpp \
	CopyMesh.cpp \
	MoveMesh.cpp \
	ExtrudeMesh.cpp \
	MergeMesh.cpp \
	Transform.cpp \
	vec_utils.cpp \
	vec_utils.hpp \
	ProjectShell.cpp \
	IntersectMesh.cpp \
	EdgeMesher.cpp \
	OneToOneSwept.cpp \
	CESets.cpp \
	CopyUtils.cpp \
	CopyUtils.hpp

# The list of header files which are to be installed
libMKalgs_la_include_HEADERS = \
	CopyGeom.hpp \
	CopyMesh.hpp \
	MoveMesh.hpp \
	ExtrudeMesh.hpp \
	MergeMesh.hpp \
	Transform.hpp \
	ProjectShell.hpp \
	IntersectMesh.hpp \
	EdgeMesher.hpp \
	OneToOneSwept.hpp \
	CESets.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils


