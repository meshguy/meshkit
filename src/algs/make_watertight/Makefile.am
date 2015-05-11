
# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libMakeWatertight.la 

# The directory in which to install headers
libMakeWatertight_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMakeWatertight_la_SOURCES = \
	arc.cpp \
	cleanup.cpp \
	gen.cpp \
        zip.cpp \
	cw_func.cpp \
        mw_func.cpp

# The list of header files which are to be installed
libMakeWatertight_la_include_HEADERS = \
        meshkit/arc.hpp \
        meshkit/cleanup.hpp \
        meshkit/gen.hpp \
        meshkit/zip.hpp \
        meshkit/cw_func.hpp \
        meshkit/mw_func.hpp

bin_PROGRAMS = make_watertight check_watertight
make_watertight_SOURCES = make_watertight.cpp 
make_watertight_CXXFLAGS = -DTEST $(CXXFLAGS)
make_watertight_LDADD = libMakeWatertight.la $(IMESH_LIBS) $(IGEOM_LIBS)

check_watertight_SOURCES = check_watertight.cpp
check_watertight_CXXFLAGS = -DTEST $(CXXFLAGS)
check_watertight_LDADD = libMakeWatertight.la $(IMESH_LIBS) $(IGEOM_LIBS)


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
