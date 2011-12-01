# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libSweep.la
         
SUBDIRS = 

# The directory in which to install headers
libSweep_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libSweep_la_SOURCES = OneToOneSwept.cpp \
        TFIMapping.cpp \
        EquipotentialSmooth.cpp \
        EquipotentialSmooth.hpp
#        HarmonicMap.cpp \
#        HarmonicMap.hpp \
#        MatrixOp.hpp \
#        MatrixOp.cpp \
#        SmoothAngleBased.hpp \
#        SmoothAngleBased.cpp

if HAVE_MESQUITE
libSweep_la_SOURCES += SweepWrapper.hpp \
        SweepWrapper.cpp \
        MeshImprove.cpp \
        SmartLaplaceWrapper.hpp \
        SmartLaplaceWrapper.cpp
endif


# The list of header files which are to be installed
nobase_libSweep_la_include_HEADERS = \
        meshkit/OneToOneSwept.hpp \
        meshkit/TFIMapping.hpp

if HAVE_MESQUITE
nobase_libSweep_la_include_HEADERS += \
        meshkit/MeshImprove.hpp 
endif

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(srcdir)/../../core \
              -I$(srcdir)/../../lemon \
              -I$(srcdir)/../../utils \
              -I$(srcdir)/../../algs \
              -I$(top_builddir) \
              -I$(top_srcdir)/utils 

if HAVE_MESQUITE
AM_CPPFLAGS += $(MESQUITE_CPPFLAGS) -DHAVE_MESQUITE
LIBS += $(MESQUITE_LIBS)

AM_CPPFLAGS += $(MESQUITE_DEFINES) \
              $(MESQUITE_CPPFLAGS)

AM_LDFLAGS += $(MESQUITE_LDFLAGS)
endif

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif
