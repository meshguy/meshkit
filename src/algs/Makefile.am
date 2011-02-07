# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitalgs.la
         
SUBDIRS = 

# The directory in which to install headers
libMeshKitalgs_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitalgs_la_SOURCES = \
	EdgeMesher.cpp \
	VertexMesher.cpp \
	SCDMesh.cpp \
	EBMesher.cpp \
	CopyMesh.cpp \
	ProjectShell.cpp
	utils/CESets.cpp \
	utils/CESets.hpp \
	utils/Transform.cpp \
	utils/Transform.hpp \
	utils/CopyUtils.cpp \
	utils/CopyUtils.hpp \
	utils/vec_utils.cpp \
	utils/vec_utils.hpp 

# The list of header files which are to be installed
nobase_libMeshKitalgs_la_include_HEADERS = \
	meshkit/EdgeMesher.hpp \
	meshkit/VertexMesher.hpp \
	meshkit/SCDMesh.hpp \
	meshkit/EBMesher.hpp \
	meshkit/CopyMesh.hpp \
	meshkit/ProjectShell.hpp


# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../core -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/ -I$(top_srcdir)/utils \
              -I$(srcdir)/../lemon

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
