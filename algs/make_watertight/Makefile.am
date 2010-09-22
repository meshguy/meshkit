@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMakeWatertight.la 

# The directory in which to install headers
libMakeWatertight_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMakeWatertight_la_SOURCES = \
	arc.cpp \
	arc.hpp \
	cleanup.cpp \
	cleanup.hpp \
	gen.cpp \
	gen.hpp \
	zip.cpp \
	zip.hpp

# The list of header files which are to be installed
libMakeWatertight_la_include_HEADERS = 

# Tests and such
check_PROGRAMS = 

TESTS = $(check_PROGRAMS)

bin_PROGRAMS = make_watertight
make_watertight_SOURCES = make_watertight.cpp 
make_watertight_CXXFLAGS = -DTEST $(CXXFLAGS)
make_watertight_LDADD = libMakeWatertight.la $(IMESH_LIBS) $(IGEOM_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils

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

libMakeWatertight_la_LIBADD = $(IMESH_LIBDIR)/libdagmc.la 
libMakeWatertight_la_DEPENDENCIES = 

dumtest:
	cat $(IMESH_LIBDIR)/libdagmc.la
	echo $(IMESH_INCLUDES)

