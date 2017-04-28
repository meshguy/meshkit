set(AssyMesher_srcs
  AssyMesher.cpp )

set(AssyMesher_headers
  meshkit/AssyMesher.hpp)

add_definitions("-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

if (MOAB_FOUND)
  add_definitions(-DMOAB)
endif ()

add_library(AssyMesher
  ${AssyMesher_srcs}
  ${AssyMesher_headers})

link_libraries(AssyMesher
  MeshKitcore)

include_directories(
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/utils
    ${CMAKE_SOURCE_DIR}/src/lemon
    ${CMAKE_SOURCE_DIR}/src/algs/Sweep
    
    )

install(
  TARGETS   "AssyMesher"
  EXPORT    meshkit
  RUNTIME   DESTINATION bin
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Runtime)
