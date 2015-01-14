set(CoreGen_srcs
  CoreGen.cpp
  copymove_parallel.cpp)

set(CoreGen_headers
  meshkit/CoreGen.hpp)

set(extra_libs)
add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")
if (MOAB_FOUND)
  add_definitions(-DMOAB)
  list(APPEND extra_libs
    iMesh)
endif ()

add_library(CoreGen STATIC
  ${CoreGen_srcs}
  ${CoreGen_headers})
target_link_libraries(CoreGen
  MeshKitcore
  ${extra_libs})
target_include_directories(CoreGen
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}")
