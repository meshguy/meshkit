# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libAlgsQslim.la

# The directory in which to install headers
libAlgsQslim_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libAlgsQslim_la_SOURCES = QslimOptions.cpp QslimDecimation.cpp \
              3D.cpp  quadrics.cpp \
	     Heap.cpp Mat4.cpp primitives.cpp

# The list of header files which are to be installed
libAlgsQslim_la_include_HEADERS = \
	QslimDecimation.h \
	QslimOptions.h 


# Tests and such
check_PROGRAMS = Qslim_test 

TESTS = $(check_PROGRAMS)

Qslim_test_SOURCES = Qslim_test.cpp 
Qslim_test_DEPENDENCIES = QslimDecimation.h libAlgsQslim.la
Qslim_test_CXXFLAGS = -DTEST $(CXXFLAGS)
Qslim_test_LDADD = libAlgsQslim.la $(IMESH_LIBS) $(IGEOM_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/utils

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


