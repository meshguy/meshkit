# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

INPUTDIR = $(top_srcdir)/rgg

# Subdirectories to build
SUBDIRS = 

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
#TESTS = $(check_PROGRAMS)
TESTS = pin217

assygen_SOURCES = assygen.cpp
assygen_DEPENDENCIES = libMKrgg.la  
assygen_CXXFLAGS = -DTEST $(CXXFLAGS)
assygen_LDADD = libMKrgg.la ../algs/libMKalgs.la ../utils/libMKutils.la $(IGEOM_LIBS)  $(IMESH_LIBS)  $(IREL_LIBS)

coregen_SOURCES = coregen.cpp 
coregen_DEPENDENCIES = ../algs/CopyMesh.hpp ../algs/MergeMesh.hpp libMKrgg.la ../algs/libMKalgs.la 
coregen_CXXFLAGS = -DTEST $(CXXFLAGS)
coregen_LDADD = libMKrgg.la ../algs/libMKalgs.la ../utils/libMKutils.la $(IGEOM_LIBS) $(IMESH_LIBS)  $(IREL_LIBS)

pin217_SOURCES = pin217.cpp utils.cpp
pin217_DEPENDENCIES = ../algs/CopyMesh.hpp ../algs/MergeMesh.hpp ../algs/libMKalgs.la 
pin217_CXXFLAGS = -DTEST $(CXXFLAGS)
pin217_LDADD = ../algs/libMKalgs.la ../utils/libMKutils.la  $(IGEOM_LIBS)  $(IMESH_LIBS)  $(IREL_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
               -I$(top_srcdir) \
               -I$(top_builddir) \
               -I$(top_srcdir)/algs \
               -I$(top_srcdir)/utils

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

