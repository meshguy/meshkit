@IGEOM_CONFIG_OPTIONS_INCLUDE@
@IMESH_CONFIG_OPTIONS_INCLUDE@
@IREL_CONFIG_OPTIONS_INCLUDE@

# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign
# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

# Subdirectories to build
SUBDIRS = . test
if ENABLE_utils
  SUBDIRS += utils
  AM_CPPFLAGS += -I$(top_srcdir)/utils
endif
if ENABLE_algs
  SUBDIRS += algs
  AM_CPPFLAGS += -I$(top_srcdir)/algs
endif
if ENABLE_rgg
  SUBDIRS += rgg
  AM_CPPFLAGS += -I$(top_srcdir)/rgg
endif
if ENABLE_src
  SUBDIRS += src
  AM_CPPFLAGS += -I$(top_srcdir)/src 
endif
if HAVE_CAMAL
  SUBDIRS += camal
  AM_CPPFLAGS += -I$(top_srcdir)/camal
endif

SUBDIRS += doc

# Some variables
AM_CPPFLAGS += -I$(top_builddir) -I$(top_srcdir) $(DEFINES)

check_PROGRAMS =  
TESTS = $(check_PROGRAMS)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = meshkit.make

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "MESHKIT_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
	echo "MESHKIT_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
