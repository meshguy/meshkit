# Things to build
lib_LTLIBRARIES = libMeshKit.la 

# Subdirectories to build
SUBDIRS = core algs

# Some variables
libMeshKit_la_LIBADD = core/libMeshKitcore.la algs/libMeshKitalgs.la
AM_CPPFLAGS += -I$(srcdir)/core -I(srcdir)/algs

# The directory in which to install headers
libMeshKit_la_includedir = $(includedir)


# The list of source files, and any header files that do not need to be installed
libMeshKit_la_SOURCES = 

# The list of header files which are to be installed
nobase_libMeshKit_la_include_HEADERS = 
