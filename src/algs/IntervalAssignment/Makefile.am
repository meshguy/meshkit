# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libIntervalAssignment.la
         
SUBDIRS = 

if HAVE_IPOPT
  AM_CPPFLAGS += $(IPOPT_CPPFLAGS)
  LDFLAGS += $(IPOPT_LDFLAGS)
endif

# The directory in which to install headers
libIntervalAssignment_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libIntervalAssignment_la_SOURCES = \
IAData.hpp \
IAInterface.cpp \
IANlp.cpp \
IANlp.hpp \
IARoundingNlp.cpp \
IARoundingNlp.hpp \
IASolution.hpp \
IASolver.cpp \
IASolver.hpp \
IASolverInt.cpp \
IASolverInt.hpp \
IASolverRelaxed.cpp \
IASolverRelaxed.hpp \
IPData.cpp \
IPData.hpp 

# The list of header files which are to be installed
nobase_libIntervalAssignment_la_include_HEADERS = \
        meshkit/IAInterface.hpp
		meshkit/IAVariable.hpp

AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core \
              -I$(top_srcdir)/src/utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/src/lemon \
              -I$(top_srcdir)/src/algs/QuadMesher \
              -I$(top_srcdir)/src/algs/Sweep \
	      -I$(top_srcdir)/src/algs/AssyGen

TESTS = test_intassign
check_PROGRAMS = $(TESTS)

LDADD = $(top_builddir)/src/libMeshKit.la

test_intassign_SOURCES = ia_main.cpp
test_intassign_DEPENDENCIES = $(LDADD)
