# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libIntervalAssignment.la
         
SUBDIRS = 

# The directory in which to install headers
libIntervalAssignment_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libIntervalAssignment_la_SOURCES = \
IAData.hpp \
IAIntWaveNlp.cpp \
IAIntWaveNlp.hpp \
IAInterface.cpp \
IAMINlp.cpp \
IAMINlp.hpp \
IANlp.cpp \
IANlp.hpp \
IASolution.hpp \
IASolver.cpp \
IASolver.hpp \
IASolverEven.cpp \
IASolverEven.hpp \
IASolverInt.cpp \
IASolverInt.hpp \
IASolverRelaxed.cpp \
IASolverRelaxed.hpp \
IPData.cpp \
IPData.hpp 

# The list of header files which are to be installed
nobase_libIntervalAssignment_la_include_HEADERS = \
        meshkit/IAInterface.hpp
		meshkit/IAVariable.hpp
