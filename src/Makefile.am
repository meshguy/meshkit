# Things to build
lib_LTLIBRARIES = libMeshKit.la 

# Subdirectories to build
SUBDIRS = core algs lemon

# Some variables
libMeshKit_la_LIBADD = core/libMeshKitcore.la algs/libMeshKitalgs.la lemon/libemon.la
AM_CPPFLAGS += -I$(srcdir)/core -I$(srcdir)/algs -I$(srcdir)/utils -I$(srcdir)/lemon -I$(srcdir)/lemon/concepts -I$(srcdir)/lemon/bits

# The directory in which to install headers
libMeshKit_la_includedir = $(includedir)

if HAVE_IGEOM
AM_CPPFLAGS += $(IGEOM_INCLUDES) $(IGEOM_CPPFLAGS)
LIBS += $(IGEOM_LIBS)
endif

if HAVE_IMESH
AM_CPPFLAGS += $(IMESH_INCLUDES) $(IMESH_CPPFLAGS)
LIBS += $(IMESH_LIBS)
endif

if HAVE_IREL
AM_CPPFLAGS += $(IREL_INCLUDES) $(IREL_CPPFLAGS)
LIBS += $(IREL_LIBS)
endif

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

# The list of source files, and any header files that do not need to be installed
libMeshKit_la_SOURCES = 

# The list of header files which are to be installed
nobase_libMeshKit_la_include_HEADERS = 
