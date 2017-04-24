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

add_library(CoreGen 
  ${CoreGen_srcs}
  ${CoreGen_headers})
link_libraries(CoreGen
  MeshKitcore
  ${extra_libs})
include_directories(CoreGen
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

install(
  TARGETS   "CoreGen"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
