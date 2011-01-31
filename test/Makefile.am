SUBDIRS = core algs

build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

DATADIR = $(top_srcdir)/data

AM_CPPFLAGS += -DDATADIR=$(DATADIR)

