# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

@IGEOM_CONFIG_OPTIONS@

# Things to build
lib_LTLIBRARIES = libMKrgg.la

# The directory in which to install headers
libMKrgg_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libMKrgg_la_SOURCES = \
	nrgen.cpp \
	io.cpp\
	clock.cpp \
	parser.cpp \
	pincell.cpp  \
	cylinder.cpp 
	

# The list of header files which are to be installed
libMKrgg_la_include_HEADERS = \
	nrgen.hpp \
	clock.hpp \
	pincell.hpp \
	cylinder.hpp \
	vectortemplate.hpp \
        matrixtemplate.hpp \
	simplearray.hpp

# Tests and such
check_PROGRAMS = assygen
TESTS = $(check_PROGRAMS)

assygen_SOURCES = assygen.cpp
assygen_DEPENDENCIES = nrgen.hpp libMKrgg.la
assygen_CXXFLAGS = -DTEST $(CXXFLAGS)
assygen_LDADD = libMKrgg.la $(IGEOM_LIBS) 

# Boilerplate stuff that doesn't depend on what the targets are
DEFS = $(DEFINES) -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)

if HAVE_IGEOM
INCLUDES += $(IGEOM_INCLUDES)
CPPFLAGS += $(IGEOM_CPPFLAGS)
endif


