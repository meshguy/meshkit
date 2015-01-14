cmake_minimum_required(VERSION 2.8.12)
project(meshkit C CXX)

option(WITH_MPI "Enable MPI support"      OFF)
option(WITH_DEBUG "Enable debug support"  OFF)
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

if (WITH_DEBUG)
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

if (WITH_MPI)
  add_definitions(-DUSE_MPI)
  find_package(MPI REQUIRED)
endif ()

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

find_package(MOAB REQUIRED)
add_definitions(-DHAVE_MOAB)
if (MOAB_HAS_FBIGEOM)
  add_definitions(-DENABLE_FBiGeom)
endif ()
if (MOAB_HAS_IMESH)
  add_definitions(-DENABLE_iMesh)
else ()
  message(FATAL_ERROR "MOAB must be built with iMesh enabled")
endif ()

find_package(CGM)
if (CGM_FOUND)
  add_definitions(-DHAVE_CGM)
  add_definitions(-DHAVE_iGeom)
  if (CGM_HAS_PARALLEL)
    add_definitions(-DHAVE_PARALLEL_CGM)
  endif ()
endif ()

find_package(Lasso)
if (Lasso_FOUND)
  add_definitions(-DHAVE_iRel)
endif ()

include(CMakeDependentOption)
cmake_dependent_option(ENABLE_UTILS "Enable utilities" ON
  "CGM_FOUND;MOAB_FOUND;MOAB_HAS_IMESH" OFF)
cmake_dependent_option(ENABLE_ALGS "Enable algorithms" ON
  "ENABLE_UTILS;Lasso_FOUND" OFF)
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
  add_subdirectory(quadmesher)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/quadmesher")
endif ()

if (ENABLE_SRC)
  add_subdirectory(src)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
endif ()

if (ENABLE_INTASSIGN)
  find_package(ipopt)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/algs/IntervalAssignment")
endif ()

include_directories(
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}")

#add_subdirectory(test)
#add_subdirectory(doc)

get_property(meshkit_targets GLOBAL
  PROPERTY meshkit_export_targets)

export(
  TARGETS ${meshkit_targets}
  FILE    "${CMAKE_CURRENT_BINARY_DIR}/MeshKitTargets.cmake")

if (WIN32)
  set(cmakedir cmake)
else ()
  set(cmakedir lib/cmake/MeshKit)
endif ()

install(
  EXPORT      meshkit
  DESTINATION "${cmakedir}"
  FILE        MeshKitTargets.cmake
  COMPONENT   development)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/MeshKitConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/MeshKitConfig.cmake"
  @ONLY)
install(
  FILES       "${CMAKE_CURRENT_BINARY_DIR}/MeshKitConfig.cmake"
  DESTINATION "${cmakedir}"
  COMPONENT   development)
