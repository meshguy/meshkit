# Things to build
lib_LTLIBRARIES = libMeshKitextern.la 

# Subdirectories to build
SUBDIRS = 

if HAVE_CAMAL
SUBDIRS += CAMAL
libMeshKitextern_la_LIBADD = CAMAL/libMeshKitCAMAL.la
AM_CPPFLAGS += -I$(srcdir)/CAMAL $(CAMAL_CPPFLAGS)
LIBS += $(CAMAL_LIBS)
endif

# The directory in which to install headers
libMeshKitextern_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKitextern_la_SOURCES = 

# The list of header files which are to be installed
nobase_libMeshKitextern_la_include_HEADERS = 
