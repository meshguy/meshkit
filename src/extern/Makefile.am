# Things to build
noinst_LTLIBRARIES = libMeshKitextern.la 

# Subdirectories to build
SUBDIRS = 

AM_CPPFLAGS += -I$(srcdir)/../core \
               -I$(srcdir)/../utils \
               -I$(srcdir)/../lemon

LIBADD =

if HAVE_CAMAL
SUBDIRS += CAMAL
LIBADD += CAMAL/libMeshKitCAMAL.la
AM_CPPFLAGS += -I$(srcdir)/CAMAL $(CAMAL_CPPFLAGS) -DHAVE_CAMAL
LIBS += $(CAMAL_LIBS)
endif

if HAVE_NETGEN
SUBDIRS += netgen
LIBADD += netgen/libMeshKitnetgen.la
AM_CPPFLAGS += -I$(srcdir)/netgen $(NETGEN_CPPFLAGS) -DHAVE_NETGEN
LIBS += $(NETGEN_LIBS)
endif

if HAVE_TETGEN
SUBDIRS += tetgen
LIBADD += tetgen/libMeshKittetgen.la
AM_CPPFLAGS += -I$(srcdir)/tetgen $(TETGEN_CPPFLAGS) -DHAVE_TETGEN
LIBS += $(TETGEN_LIBS)
endif

if HAVE_TRIANGLE
SUBDIRS += triangle
LIBADD += triangle/libMeshKittriangle.la
AM_CPPFLAGS += -I$(srcdir)/triangle $(TRIANGLE_CPPFLAGS) -DHAVE_TRIANGLE
LIBS += $(TRIANGLE_LIBS)
endif

if HAVE_MESQUITE
SUBDIRS += mesquite
LIBADD += mesquite/libMeshKitmesquite.la
AM_CPPFLAGS += -I$(srcdir)/mesquite $(MESQUITE_CPPFLAGS)  -DHAVE_MESQUITE
LIBS += $(MESQUITE_LIBS)
endif

# The directory in which to install headers
libMeshKitextern_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitextern_la_SOURCES = register_extern.cpp

# The list of header files which are to be installed
nobase_libMeshKitextern_la_include_HEADERS = 

libMeshKitextern_la_LIBADD = $(LIBADD)

