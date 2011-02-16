# Things to build
lib_LTLIBRARIES = libMeshKit.la 

# Subdirectories to build
SUBDIRS = core algs utils lemon extern

# Some variables
libMeshKit_la_LIBADD = core/libMeshKitcore.la \
                       algs/libMeshKitalgs.la \
                       lemon/libemon.la \
                       utils/libMeshKitutils.la \
                       extern/libMeshKitextern.la
AM_CPPFLAGS += -I$(srcdir)/core -I$(srcdir)/algs -I$(srcdir)/utils -I$(srcdir)/lemon -I$(srcdir)/lemon/concepts -I$(srcdir)/lemon/bits

# The directory in which to install headers
libMeshKit_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMeshKit_la_SOURCES = 

# The list of header files which are to be installed
nobase_libMeshKit_la_include_HEADERS = 
