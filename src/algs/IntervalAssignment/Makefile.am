# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libIntervalAssignment.la
         
SUBDIRS = 

if HAVE_INTASSIGN
  AM_CPPFLAGS += $(IPOPT_CPPFLAGS)
  LDFLAGS += $(IPOPT_LDFLAGS)
endif

# The directory in which to install headers
libIntervalAssignment_la_includedir = $(includedir)

libIntervalAssignment_la_LIBADD = $(IPOPT_LIBS)

# The list of source files, and any header files that do not need to be installed
libIntervalAssignment_la_SOURCES = \
IADataBuilder.cpp \
IAInterface.cpp \
IAIntWaveNlp.cpp \
IAIntParabolaNlp.cpp \
IANlp.cpp \
IARoundingNlp.cpp \
IASolver.cpp \
IASolverTool.cpp \
IASolverToolInt.cpp \
IASolverInt.cpp \
IASolverRelaxed.cpp \
IAWeights.cpp \
IPData.cpp \
IPBend.cpp \
IASolverBend.cpp \
IABendNlp.cpp

# The list of header files which are to be installed
nobase_libIntervalAssignment_la_include_HEADERS = \
meshkit/IABendNlp.hpp \
meshkit/IADataBuilder.hpp \
meshkit/IAData.hpp \
meshkit/IAIntCosNlp.hpp \
meshkit/IAInterface.hpp \
meshkit/IAIntParabolaNlp.hpp \
meshkit/IAIntWaveNlp.hpp \
meshkit/IAMilp.hpp \
meshkit/IAMINlp.hpp \
meshkit/IANlp.hpp \
meshkit/IARoundingFar3StepNlp.hpp \
meshkit/IARoundingFarNlp.hpp \
meshkit/IARoundingHeuristicMINLP.hpp \
meshkit/IARoundingNlp.hpp \
meshkit/IASolution.hpp \
meshkit/IASolverBend.hpp \
meshkit/IASolverEven.hpp \
meshkit/IASolver.hpp \
meshkit/IASolverInt.hpp \
meshkit/IASolverRelaxed.hpp \
meshkit/IASolverTool.hpp \
meshkit/IASolverToolInt.hpp \
meshkit/IAVariable.hpp \
meshkit/IAWeights.hpp \
meshkit/IPBend.hpp \
meshkit/IPData.hpp


AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core \
              -I$(top_srcdir)/src/utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/src/lemon \
              -I$(top_srcdir)/src/algs/QuadMesher \
              -I$(top_srcdir)/src/algs/Sweep \
	      -I$(top_srcdir)/src/algs/AssyGen

LDADD = $(top_builddir)/src/libMeshKit.la

