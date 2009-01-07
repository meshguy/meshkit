# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign
# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

# Subdirectories to build
SUBDIRS = 
if ENABLE_algs
  SUBDIRS += algs
  INCLUDES += -I$(top_builddir)/algs
endif
if ENABLE_xforms
  SUBDIRS += xforms
  INCLUDES += -I$(top_builddir)/xforms
endif
if ENABLE_utils
  SUBDIRS += utils
  INCLUDES += -I$(top_builddir)/utils
endif

# Some variables
DEFS = $(DEFINES) 
INCLUDES += -I$(top_builddir) 

check_PROGRAMS = pin217 sixth
TESTS = $(check_PROGRAMS)

pin217_SOURCES = pin217.cpp
pin217_INCLUDES = algs utils
pin217_DEPENDENCIES = algs/CopyMesh.hpp algs/MergeMesh.hpp algs/libMKalgs.la utils/libMKutils.la
pin217_CXXFLAGS = $(CXXFLAGS)
pin217_LDADD = algs/libMKalgs.la utils/libMKutils.la $(IMESH_LIBS)

sixth_SOURCES = sixth.cpp
sixth_INCLUDES = algs utils
sixth_DEPENDENCIES = algs/CopyMesh.hpp algs/MergeMesh.hpp algs/libMKalgs.la utils/libMKutils.la
sixth_CXXFLAGS = $(CXXFLAGS)
sixth_LDADD = algs/libMKalgs.la utils/libMKutils.la $(IMESH_LIBS)

if HAVE_IGEOM
INCLUDES += $(iGeom_INCLUDES)
CPPFLAGS += $(iGeom_CPPFLAGS)
LIBS += $(IGEOM_LIBS)
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
LIBS += $(IMESH_LIBS)
#CXXFLAGS += $(IMESH_CXXFLAGS)
#CFLAGS += $(IMESH_CFLAGS)
#LDFLAGS += $(IMESH_LDFLAGS)
endif

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = meshkit.make

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "MESHKIT_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
	echo "MESHKIT_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
