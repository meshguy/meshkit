
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
	zip.hpp \
	cw_func.cpp \
	cw_func.hpp \
	mw_func.cpp \
	mw_func.hpp

# The list of header files which are to be installed
libMakeWatertight_la_include_HEADERS = 

# Tests and such
TESTS = 

bin_PROGRAMS = make_watertight
make_watertight_SOURCES = make_watertight.cpp 
make_watertight_CXXFLAGS = -DTEST $(CXXFLAGS)
make_watertight_LDADD = libMakeWatertight.la $(IMESH_LIBS) $(IGEOM_LIBS)


libMakeWatertight_la_LIBADD = -ldagmc
libMakeWatertight_la_DEPENDENCIES = 

dumtest:
	cat $(IMESH_LIBDIR)/libdagmc.la
	echo $(IMESH_INCLUDES)

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
               -I$(top_srcdir) \
               -I$(top_builddir) \
               -I$(top_srcdir)/algs \
               -I$(top_srcdir)/utils

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

LDADD = $(top_builddir)/src/libMeshKit.la
