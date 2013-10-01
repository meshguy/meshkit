# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libQslim.la

# The directory in which to install headers
libQslim_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libQslim_la_SOURCES = \
	3D.cpp    Heap.hpp         Mat4.h         \
	QslimDecimation.hpp  std.h 3D.h   QslimMesher.cpp   \
	Vec3.h Array.h   QslimOptions.cpp     Vec4.h \
	Buffer.h    primitives.cpp     defs.h    \
	primitives.h         quadrics.cpp   Heap.cpp  Mat4.cpp\
	QslimDecimation.cpp  quadrics.h

# The list of header files which are to be installed (from algs/meshkit)
nobase_libQslim_la_include_HEADERS = \
	meshkit/QslimMesher.hpp \
	meshkit/QslimOptions.hpp 

# Boilerplate stuff that doesn't depend on what the targets are
AM_CPPFLAGS += $(DEFINES) -DSRCDIR=$(srcdir) \
              -I$(top_srcdir)/src/core \
              -I$(top_srcdir)/src/utils \
              -I$(top_srcdir) \
              -I$(top_builddir) \
              -I$(top_srcdir)/src/lemon

LDADD = $(top_builddir)/src/libMeshKit.la
