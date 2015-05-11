# Subdirectories to build
SUBDIRS = 

if HAVE_CAMAL
SUBDIRS += CAMAL
endif

# CAMAL dependency should go once we have another trimesher.
if HAVE_NETGEN
if HAVE_CAMAL
SUBDIRS += netgen
endif
endif

if HAVE_MESQUITE
SUBDIRS += mesquite
endif

if HAVE_TRIANGLE
SUBDIRS += triangle
endif
