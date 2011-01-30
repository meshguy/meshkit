SUBDIRS = core algs

build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

if HAVE_MOAB
AM_CPPFLAGS += -DMOAB
endif

