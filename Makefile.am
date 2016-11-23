
# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign
# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I config

# Subdirectories to build
SUBDIRS = . python
if ENABLE_quadmesher
  AM_CPPFLAGS += -I$(top_srcdir)/src/algs/QuadMesher
endif
if ENABLE_src
  SUBDIRS += src
  AM_CPPFLAGS += -I$(top_srcdir)/src 
endif
if ENABLE_intassign
  AM_CPPFLAGS += -I$(top_srcdir)/src/algs/IntervalAssignment $(IPOPT_CPPFLAGS)
  AM_LDFLAGS += $(IPOPT_LDFLAGS)
  LIBS += $(IPOPT_LIBS)
endif

SUBDIRS += test doc tools

# Some variables
AM_CPPFLAGS += -I$(top_builddir) -I$(top_srcdir) $(DEFINES)

# check_PROGRAMS =
# TESTS = $(check_PROGRAMS)

EXTRA_DIST = MKDefines.h data tools/rgg_test_dirs test/algs/postbl_examples

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = meshkit.make

toolsdir = $(prefix)/examples/
tools_DATA = common.mk

# By default, meshkit.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "MESHKIT_LIBS=-L${libdir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
	echo "MESHKIT_INCLUDES=-I${includedir}" >> $(DESTDIR)$(cfgdir)/meshkit.make
	chmod +rw $(prefix)	
	mkdir -p $(prefix)/examples/rgg/
	cp -r $(top_srcdir)/tools/rgg_test_dirs/* $(prefix)/examples/rgg/
	mkdir -p $(prefix)/examples/postbl/
	cp -r $(top_srcdir)/test/algs/postbl_examples/* $(prefix)/examples/postbl/
	
uninstall-local:
	chmod -R +rw $(prefix)/examples/
	rm -rf $(prefix)/examples/
	
	
	
