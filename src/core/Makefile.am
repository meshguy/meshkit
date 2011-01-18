# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

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
	ModelEnt.cpp \
	SizingFunction.cpp

# The list of header files which are to be installed
nobase_libMeshKitcore_la_include_HEADERS = \
	meshkit/Error.hpp \
	meshkit/MeshOp.hpp \
	meshkit/MeshOpFactory.hpp \
	meshkit/MeshScheme.hpp \
	meshkit/ModelEnt.hpp \
	meshkit/SizingFunction.hpp

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS) -DCGM
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES) $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
AM_CPPFLAGS += $(IREL_INCLUDES) $(IREL_CPPFLAGS)
endif

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
