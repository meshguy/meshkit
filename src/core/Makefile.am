# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMeshKitcore.la
         
# Some variables
AM_CPPFLAGS += -DSRCDIR=$(srcdir)

# The directory in which to install headers
libMeshKitcore_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitcore_la_SOURCES = \
	MeshOp.cpp \
	MeshOpFactory.cpp \
	MeshScheme.cpp \
	MKCore.cpp \
	MKGraph.cpp \
	ModelEnt.cpp \
	SizingFunction.cpp \
	VecUtil.cpp

# The list of header files which are to be installed
nobase_libMeshKitcore_la_include_HEADERS = \
	meshkit/Error.hpp \
	meshkit/GraphNode.hpp \
	meshkit/MeshOp.hpp \
	meshkit/MeshOpFactory.hpp \
	meshkit/MeshScheme.hpp \
	meshkit/MKCore.hpp \
	meshkit/MKGraph.hpp \
	meshkit/ModelEnt.hpp \
	meshkit/SizingFunction.hpp \
	meshkit/VecUtil.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils \
              -I$(top_srcdir)/src/lemon

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
