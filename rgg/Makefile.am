# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKrgg.la

# The directory in which to install headers
libMKrgg_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKrgg_la_SOURCES = \
	nrgen.cpp \
	io.cpp \
	clock.cpp \
	parser.cpp \
	pincell.cpp  \
	cylinder.cpp \
	crgen.cpp \
	copymove.cpp \
	utils.cpp

# The list of header files which are to be installed
libMKrgg_la_include_HEADERS = \
	nrgen.hpp \
	clock.hpp \
	pincell.hpp \
	cylinder.hpp \
	vectortemplate.hpp \
	matrixtemplate.hpp \
	crgen.hpp \
	utils.hpp 

# Tests and such
check_PROGRAMS = assygen coregen pin217 
TESTS = $(check_PROGRAMS)

assygen_SOURCES = assygen.cpp
assygen_DEPENDENCIES = libMKrgg.la
assygen_CXXFLAGS = -DTEST $(CXXFLAGS)
assygen_LDADD = libMKrgg.la $(IGEOM_LIBS) 

coregen_SOURCES = coregen.cpp 
coregen_DEPENDENCIES = ../algs/CopyMesh.hpp ../algs/MergeMesh.hpp libMKrgg.la ../algs/libMKalgs.la 
coregen_CXXFLAGS = -DTEST $(CXXFLAGS)
coregen_LDADD = libMKrgg.la ../algs/libMKalgs.la ../utils/libMKutils.la $(IGEOM_LIBS) $(IMESH_LIBS)

pin217_SOURCES = pin217.cpp utils.cpp
pin217_DEPENDENCIES = ../algs/CopyMesh.hpp ../algs/MergeMesh.hpp ../algs/libMKalgs.la 
pin217_CXXFLAGS = -DTEST $(CXXFLAGS)
pin217_LDADD = ../algs/libMKalgs.la ../utils/libMKutils.la $(IMESH_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/algs -I$(top_srcdir)/utils

if HAVE_IGEOM
INCLUDES += $(IGEOM_INCLUDES)
CPPFLAGS += $(IGEOM_CPPFLAGS)
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
endif

if HAVE_MOAB
CPPFLAGS += -DMOAB
endif

