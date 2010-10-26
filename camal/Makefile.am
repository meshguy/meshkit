# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@IREL_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKcamal.la libSmooth.la

# The directory in which to install headers
libMKcamal_la_includedir = $(includedir)
libSmooth_la_includedir = $(includedir)

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

# The list of source files, and any header files that do not need to be installed
libSmooth_la_SOURCES = \
	SmoothFaceEval.cpp \
	CamalPaveDriver.cpp \
	CAMALSizeEval.cpp \
	SmoothCurveEval.cpp \
	SmoothVertex.cpp  \
	SmoothBase.cpp

# The list of header files which are to be installed
libSmooth_la_include_HEADERS = \
	SmoothFaceEval.hpp \
	CamalPaveDriver.hpp \
	CAMALSizeEval.hpp \
	SmoothCurveEval.hpp \
	SmoothVertex.hpp \
	SmoothBase.hpp
	

# Tests
check_PROGRAMS = Camal_test IceSheet_test SmoothIceSheet_test addGeoSets Create2FaceShell \
					NewSmooth_test ComputeThickness \
                                         extrudeQuads

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

SmoothIceSheet_test_SOURCES = smooth_test.cpp
SmoothIceSheet_test_DEPENDENCIES = libSmooth.la
SmoothIceSheet_test_CXXFLAGS = -DTEST $(CXXFLAGS)
SmoothIceSheet_test_LDADD = libSmooth.la $(CAMAL_LIBS) $(IMESH_LIBS)

addGeoSets_SOURCES = Geometrize.cpp
addGeoSets_DEPENDENCIES = libSmooth.la
addGeoSets_CXXFLAGS = -DTEST $(CXXFLAGS)
addGeoSets_LDADD = libSmooth.la $(IMESH_LIBS) 

Create2FaceShell_SOURCES = CreateShell.cpp
Create2FaceShell_DEPENDENCIES =  
Create2FaceShell_CXXFLAGS = -DTEST $(CXXFLAGS)
Create2FaceShell_LDADD = $(IMESH_LIBS)

NewSmooth_test_SOURCES = newSmoothTest.cpp
NewSmooth_test_DEPENDENCIES = libMKcamal.la
NewSmooth_test_CXXFLAGS = -DTEST $(CXXFLAGS)
NewSmooth_test_LDADD = libMKcamal.la $(CAMAL_LIBS) $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

ComputeThickness_SOURCES = ComputeThickness.cpp
ComputeThickness_DEPENDENCIES =  
ComputeThickness_CXXFLAGS = $(CXXFLAGS)
ComputeThickness_LDADD =  $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)

extrudeQuads_SOURCES = extrudeQuads.cpp
extrudeQuads_DEPENDENCIES =  
extrudeQuads_CXXFLAGS = $(CXXFLAGS)
extrudeQuads_LDADD =  $(IMESH_LIBS) $(IGEOM_LIBS) $(IREL_LIBS)


# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS)
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES)  $(IMESH_CPPFLAGS)
endif

if HAVE_IREL
AM_CPPFLAGS += $(IREL_INCLUDES) $(IREL_CPPFLAGS)
endif

if HAVE_CAMAL
AM_CPPFLAGS += -DCAMAL_VERSION=$(CAMAL_VERSION) $(CAMAL_CPPFLAGS)
LDFLAGS += $(CAMAL_LDFLAGS)
endif
