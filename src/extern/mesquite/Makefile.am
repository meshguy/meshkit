noinst_LTLIBRARIES = libMeshKitmesquite.la
libMeshKitmesquite_la_SOURCES = MesquiteOpt.cpp
nobase_include_HEADERS = meshkit/MesquiteOpt.hpp
AM_CPPFLAGS += -I$(top_srcdir)/src/utils \
               -I$(top_srcdir)/src/core \
               -I$(top_srcdir)/src/lemon \
               $(MESQUITE_DEFINES) \
               $(MESQUITE_CPPFLAGS) 
AM_LDFLAGS += $(MESQUITE_LDFLAGS)
libMeshKitmesquite_la_LIBADD = $(MESQUITE_LIBS)

