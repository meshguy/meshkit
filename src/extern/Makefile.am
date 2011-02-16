# Things to build
noinst_LTLIBRARIES = libMeshKitextern.la 

# Subdirectories to build
SUBDIRS = 

AM_CPPFLAGS += -I$(srcdir)/../core \
               -I$(srcdir)/../utils \
               -I$(srcdir)/../lemon


if HAVE_CAMAL
SUBDIRS += CAMAL
libMeshKitextern_la_LIBADD = CAMAL/libMeshKitCAMAL.la
AM_CPPFLAGS += -I$(srcdir)/CAMAL $(CAMAL_CPPFLAGS) -DHAVE_CAMAL
LIBS += $(CAMAL_LIBS)
endif

# The directory in which to install headers
libMeshKitextern_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitextern_la_SOURCES = register_extern.cpp

# The list of header files which are to be installed
nobase_libMeshKitextern_la_include_HEADERS = 
