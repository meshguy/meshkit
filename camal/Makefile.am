# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKcamal.la

# The directory in which to install headers
libMKcamal_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKcamal_la_SOURCES = \
	camal_interface.cpp \
	camel.cpp \
	CAMALGeomEval.cpp\
	CAMALSizeEval.cpp

# The list of header files which are to be installed
libMKcamal_la_include_HEADERS = \
	camal_interface.hpp \
	camel.hpp \
	CAMALGeomEval.hpp \
	CAMALSizeEval.hpp

# Tests
check_PROGRAMS = Camal_test IceSheet_test

if HAVE_CAMAL
TESTS = $(check_PROGRAMS)
endif

Camal_test_SOURCES = main.cpp
Camal_test_DEPENDENCIES = libMKcamal.la
Camal_test_CXXFLAGS = -DTEST $(CXXFLAGS)
Camal_test_LDADD = libMKcamal.la $(CAMAL_LIBS) $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

IceSheet_test_SOURCES = icesheet_test.cpp
IceSheet_test_DEPENDENCIES = libMKcamal.la
IceSheet_test_CXXFLAGS = -DTEST $(CXXFLAGS)
IceSheet_test_LDADD = libMKcamal.la $(CAMAL_LIBS) $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
INCLUDES += $(IGEOM_INCLUDES)
CPPFLAGS += $(IGEOM_CPPFLAGS)
endif

if HAVE_IMESH
INCLUDES += $(IMESH_INCLUDES)
CPPFLAGS += $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
INCLUDES += $(IREL_INCLUDES)
CPPFLAGS += $(IREL_CPPFLAGS)
endif

if HAVE_CAMAL
CPPFLAGS += -DCAMAL_VERSION=$(CAMAL_VERSION) $(CAMAL_CPPFLAGS)
LDFLAGS += $(CAMAL_LDFLAGS)
endif