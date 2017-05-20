cmake_minimum_required( VERSION 2.8.12 )
SET(PACKAGE_NAME "MeshKit")
SET(PACKAGE_VERSION "1.4.2pre")

if (NOT WIN32)
  #This all breaks on windows.
  SET(CMAKE_Fortran_COMPILER_INIT ${CMAKE_GENERATOR_FC})
  SET(CMAKE_CXX_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_C_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_C_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_C_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_Fortran_FLAGS_INIT "-fPIC")
  SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2")
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS_INIT "")
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS_INIT "")
endif()

project( MeshKit C CXX)

include(GNUInstallDirs)

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})

option(ENABLE_MPI "Enable MPI support"      OFF)
option(ENABLE_DEBUG "Enable debug support"  OFF)
option(BUILD_SHARED_LIBS "Build shard libraries" OFF)

find_program(decompress_program
  NAMES gunzip gzip zcat
  DOC   "decompression program")
set(have_decompress OFF)
if (decompress_program)
  get_filename_component(zcat_name "${decompress_program}" NAME_WE)
  if (zcat_name STREQUAL "gunzip")
    set(ZCAT_COMMAND "${decompress_program} -c")
    set(have_decompress ON)
  elseif (zcat_name STREQUAL "gzip")
    set(ZCAT_COMMAND "${decompress_program} -cd")
    set(have_decompress ON)
  elseif (zcat_name STREQUAL "zcat")
    set(ZCAT_COMMAND "${decompress_program}")
    set(have_decompress ON)
  else ()
    message(WARNING "Unrecognized decompression executable: ${decompress_program} (${zcat_name})")
  endif ()
endif ()
if (NOT have_decompress)
  message(WARNING "Cannot run tests; no way to uncompress input files")
endif ()

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(vsnprintf "stdarg.h" HAVE_VSNPRINTF)
if (HAVE_VSNPRINTF)
  add_definitions(-DHAVE_VSNPRINTF)
endif ()

# Use APIs which allow files >2GB to be used.
add_definitions(-D_FILE_OFFSET_BITS=64)

set(VERSION_MAJOR  1)
set(VERSION_MINOR  3)
set(VERSION_PATCH  1pre)
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")
if (NOT VERSION_PATCH STREQUAL "")
  set(VERSION "${VERSION}.${VERSION_PATCH}")
else ()
  if (VERSION_MINOR STREQUAL 99)
    set(VERSION "${VERSION} (Alpha)")
  else ()
    set(VERSION "${VERSION} (Beta)")
  endif ()
endif ()

include(CheckIncludeFiles)
check_include_files("ieeefp.h" HAVE_IEEEFP_H)
if (HAVE_IEEEFP_H)
  add_definitions(-DHAVE_IEEEFP_H)
endif ()

if (ENABLE_DEBUG)
  check_include_files("valgrind/memcheck.h" HAVE_VALGRIND)
  if (HAVE_VALGRIND)
    add_definitions(-DVALGRIND)
  endif ()
endif ()

check_include_files("cstddef" HAVE_CSTDDEF)
check_include_files("stddef.h" HAVE_STDDEF_H)
check_include_files("stdio.h" HAVE_STDIO_H)
check_include_files("cstdio" HAVE_CSTDIO)

# TODO: doxygen

# check for MPI package
set (MESHKIT_HAVE_MPI OFF CACHE INTERNAL "Found necessary MPI components. Configure MeshKit with it." )
if ( ENABLE_MPI )
  set (MPI_C_COMPILER ${CC} )
  set (MPI_CXX_COMPILER ${CXX} )
  set (MPI_Fortran_COMPILER ${FC} )
  find_package( MPI REQUIRED )
  # CMake FindMPI script is sorely lacking:
  if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
    set( MPI_FOUND ON )
    add_definitions(-DUSE_MPI)
  endif ( MPI_LIBRARY AND MPI_INCLUDE_PATH )

  if ( MPI_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    set( MESHKIT_HAVE_MPI ON )
  endif ( MPI_FOUND )
endif ( ENABLE_MPI )
include (config/CheckCompilerFlags.cmake)


# TODO: tetgen stuff
# TODO: netgen stuff
# TODO: triangle stuff

if (ENABLE_INTASSIGN)
  # TODO: ipopt stuff
endif ()

# TODO: mesquite stuff
# TODO: lpsolver stuff
# TODO: armadillo stuff
# TODO: openblas stuff

find_package(CGM REQUIRED)
if (CGM_FOUND)
  add_definitions(-DHAVE_CGM)
  add_definitions(-DHAVE_iGeom)
  if (CGM_HAS_PARALLEL)
    add_definitions(-DHAVE_PARALLEL_CGM)
  endif ()
endif ()


find_file(MOAB_CMAKE_CFG MOABConfig.cmake
    HINTS ${MOAB_DIR} ${MOAB_DIR}/lib/cmake/MOAB
)

include (${MOAB_CMAKE_CFG})
#find_package(MOAB REQUIRED)
if (MOAB_FOUND)
    # Output details about CGM configuration
    message("Found MOAB: ${MOAB_DIR}")
    message("---   MOAB configuration ::")
    message("        Include Directory       : ${MOAB_INCLUDE_DIRS}")
    message("        Library Directory       : ${MOAB_LIBRARIES}")    
endif(MOAB_FOUND)
add_definitions(-DHAVE_MOAB)
if (ENABLE_FBIGEOM)
  add_definitions(-DENABLE_FBiGeom)
endif ()
if (ENABLE_IMESH)
  add_definitions(-DENABLE_iMesh)
else ()
  message(FATAL_ERROR "MOAB must be built with iMesh enabled")
endif ()

include(CMakeDependentOption)
cmake_dependent_option(ENABLE_UTILS "Enable utilities" ON
  "CGM_FOUND;MOAB_FOUND;MOAB_HAS_IMESH" OFF)
cmake_dependent_option(ENABLE_ALGS "Enable algorithms" ON
  "ENABLE_UTILS;TARGET iRel" OFF)
option(ENABLE_QUADMESHER "Enable quadmesher"  OFF)
option(ENABLE_RGG "Enable rgg"                ON)
option(ENABLE_SRC "Enable src"                ON)
option(ENABLE_PYTHON "Enable intassign"       OFF)
option(ENABLE_INTASSIGN "Enable intassign"    OFF)

if (ENABLE_QUADMESHER)
  add_definitions(-DHAVE_QUADMESHER)
endif ()

if (ENABLE_PYTHON OR TRUE) # XXX: unconditional in autotools?
  #add_subdirectory(python)
endif ()

if (ENABLE_QUADMESHER)
  add_subdirectory(src/algs/quadmesher)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/quadmesher)
endif ()

if (ENABLE_SRC)
  add_subdirectory(src)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
endif ()

if (ENABLE_INTASSIGN)
  find_package(ipopt)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/algs/IntervalAssignment)
endif ()

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

#add_subdirectory(test)
#add_subdirectory(doc)
add_subdirectory(tools)
