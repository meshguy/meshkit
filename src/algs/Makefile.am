# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMeshKitalgs.la

# The directory in which to install headers
libMeshKitalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitalgs_la_SOURCES = \
	register_algs.cpp \
	EdgeMesher.cpp \
	VertexMesher.cpp \
	SCDMesh.cpp \
	EBMesher.cpp \
	CopyMesh.cpp \
	CESets.cpp \
	Transform.cpp \
	CopyUtils.cpp \
	CopyUtils.hpp \
	ProjectShell.cpp \
	vec_utils.cpp \
	vec_utils.hpp \
	OneToOneSwept.cpp

# The list of header files which are to be installed
nobase_libMeshKitalgs_la_include_HEADERS = \
	meshkit/EdgeMesher.hpp \
	meshkit/VertexMesher.hpp \
	meshkit/SCDMesh.hpp \
	meshkit/EBMesher.hpp \
	meshkit/CopyMesh.hpp \
	meshkit/ProjectShell.hpp \
	meshkit/CESets.hpp \
	meshkit/Transform.hpp \
	meshkit/OneToOneSwept.hpp


# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core -I$(srcdir)/../utils -I$(top_srcdir) \
              -I$(top_builddir) -I$(top_srcdir)/ -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon

